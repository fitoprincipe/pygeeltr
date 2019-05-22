# -*- coding: utf-8 -*-

from geetools import tools, bitreader
import math
from collections import namedtuple

import ee


def statistics(collection, band=None, suffix='_fit', skip_outliers=True):
    """ Compute sum of squares and R2 for values of every pixel in a
    ImageCollection that has already be fitted

    :param collection: collection
    :type collection: ee.ImageCollection
    :param band: name of the band that has been fitted. If None, the first
        band will be used
    :type band: str
    :param suffix: suffix of band that holds the fitted values
    :type suffix: str
    :param skip_outliers: outliers (in difference) are not included for
        computing
    :type skip_outliers: bool
    :return:
            :ssres: residual sum of square in an ee.Image
            :r2: coefficient of determination (r2)
    :rtype: namedtuple
    """

    if not band:
        first = ee.Image(collection.first()).bandNames().get(0)
        band = ee.String(first).getInfo()

    fitted_band = band+suffix

    # Band mean
    justbands = collection.select([band, fitted_band])

    def diff(img):
        nbr = img.select([0])
        fit = img.select([1])
        dif = fit.subtract(nbr).abs().select([0], ['diff'])
        return img.addBands(dif)

    justbands = justbands.map(diff)

    # Stats
    mean = justbands.mean()
    median = justbands.median().select(['diff'], ['median'])
    stddev = justbands.reduce(ee.Reducer.stdDev()).select(['diff_stdDev'],
                                                          ['stddev'])
    interval_min = median.subtract(stddev.multiply(2)).select([0], ['min'])
    interval_max = median.add(stddev.multiply(2)).select([0], ['max'])

    # SSTOT
    sstot_image = tools.image.empty(0, ['sstot'])

    def sstot(img, sstoti):
        sstoti = ee.Image(sstoti) # cast
        mean_i = mean.select([band])
        val = img.select([band])
        diff_mean = val.subtract(mean_i)
        # excude outlier

        diff = img.select('diff')
        outlier = diff.lte(interval_min).Or(diff.gte(interval_max)).select([0], ['out'])

        # compute
        compute = diff_mean.pow(2).select([0], ['sstot'])

        # skip outliers
        if skip_outliers:
            compute = compute.updateMask(outlier.Not()).unmask()

        return sstoti.add(compute)

    sstot_image = ee.Image(justbands.iterate(sstot, sstot_image))

    # SSRES
    ssres_image = tools.image.empty(0, ['ssres'])

    def ssres(img, ssresi):
        ssresi = ee.Image(ssresi) # cast

        diff = img.select('diff')
        outlier = diff.lte(interval_min).Or(diff.gte(interval_max)).select([0], ['out'])

        compute = diff.pow(2).select([0], ['ssres'])

        # skip outliers
        if skip_outliers:
            compute = compute.updateMask(outlier.Not()).unmask()

        return ssresi.add(compute)

    ssres_image = ee.Image(justbands.iterate(ssres, ssres_image))

    # DIVIDE
    division = ssres_image.divide(sstot_image).select([0], ['r2'])

    # 1-division
    r2 = tools.image.empty(1, ['r2']).subtract(division).toFloat()

    result = namedtuple("Statistics", ["ssres", "r2"])

    return result(ssres_image, r2)


class LandTrendr(object):
    """
    LandTrendr (Landsat-based detection of Trends in Disturbance and
    Recovery) is an algorithm that temporally segments a time-series of
    images by extracting the spectral trajectories of change over time.

    **Parameters:**

    - timeseries (ImageCollection):

      Collection from which to extract trends (it's assumed thateach
      image in the collection represents one year). The first band is
      usedto find breakpoints, and all subsequent bands are fitted
      using those breakpoints.

    - fit_band (str):

      band to use for regression

    **Optionals:**

    - area (ee.Geometry): area to compute LandTrendr. If `None`, takes the
      area of the first image

    - date_measure (str): What date represent every image in the collection?.
      Defaults to `year`. Can be `month`, `day`, etc

    **Originals:**

    - maxSegments (Integer):

      Maximum number of segments to be fitted on the time series.

    - spikeThreshold (Float, default: 0.9):

      Threshold for dampening the spikes (1.0 means no dampening).

    - vertexCountOvershoot (Integer, default: 3):

      The inital model can overshoot the maxSegments + 1 vertices by this
      amount. Later, it will be prunned down to maxSegments + 1.

    - preventOneYearRecovery (Boolean, default: false):
      Prevent segments that represent one year recoveries.

    - recoveryThreshold (Float, default: 0.25):
      If a segment has a recovery rate faster than 1/recoveryThreshold
      (in years), then the segment is disallowed.

    - pvalThreshold (Float, default: 0.1):
      If the p-value of the fitted model exceeds this threshold, then the
      current model is discarded and another one is fitted using the
      Levenberg-Marquardt optimizer.

    - bestModelProportion (Float, default: 1.25):
      Takes the model with most vertices that has a p-value that is at
      most this proportion away from the model with lowest p-value.

    - minObservationsNeeded (Integer, default: 6):
      Min observations needed to perform output fitting.
    """

    def __init__(self, timeseries, fit_band, area=None,
                 date_measure='year', **kwargs):

        self.timeSeries = timeseries
        self.fit_band = fit_band
        self.date_measure = date_measure

        self.maxSegments = kwargs.get("maxSegments", 4)
        self.spikeThreshold = kwargs.get("spikeThreshold", 0.9)
        self.vertexCountOvershoot = kwargs.get("vertexCountOvershoot", 3)
        self.preventOneYearRecovery = kwargs.get("preventOneYearRecovery",
                                                 False)
        self.recoveryThreshold = kwargs.get("recoveryThreshold", 0.25)
        self.pvalThreshold = kwargs.get("pvalThreshold", 0.1)
        self.bestModelProportion = kwargs.get("bestModelProportion", 1.25)
        self.minObservationsNeeded = kwargs.get("minObservationsNeeded", 6)

        # Axes
        self.year_axis = 1
        self.band_axis = 0

        if area:
            self.area = area
        else:
            self.area = ee.Image(self.timeSeries.first()).geometry()

        self._core = None
        self._breakdown = None
        self._slope = None
        self._statistics = None
        self._date_range = None
        self._date_range_bitreader = None
        self._region = None

    @classmethod
    def Principe(cls, timeseries, index, area=None):
        """ Factory Method to create a LandTrendr class with params defined
        by Rodrigo E. Principe (fitoprincipe82@gmail.com)
        """
        newobj = cls(timeseries, index, area,
                     maxSegments=4,
                     spikeThreshold=0.1,
                     vertexCountOvershoot=0,
                     preventOneYearRecovery=False,
                     recoveryThreshold=0.9,
                     pvalThreshold=0.9,
                     bestModelProportion=0.1,
                     minObservationsNeeded=6)

        return newobj

    @classmethod
    def Kennedy(cls, timeseries, index, area=None):
        """ Factory Method to create a LandTrendr class with params defined
        by Kennedy (original)
        """
        newobj = cls(timeseries, index, area,
                     maxSegments=4,
                     spikeThreshold=0.9,
                     vertexCountOvershoot=3,
                     preventOneYearRecovery=False,
                     recoveryThreshold=0.25,
                     pvalThreshold=0.05,
                     bestModelProportion=0.75,
                     minObservationsNeeded=6)

        return newobj

    @classmethod
    def Liang(cls, timeseries, index, area=None):
        """ Factory Method to create a LandTrendr class with params defined
        by Liang
        """
        newobj = cls(timeseries, index, area,
                     maxSegments=4,
                     spikeThreshold=0.9,
                     vertexCountOvershoot=0,
                     preventOneYearRecovery=False,
                     recoveryThreshold=0.25,
                     pvalThreshold=0.1,
                     bestModelProportion=0.75,
                     minObservationsNeeded=6)

        return newobj

    @property
    def region(self):
        """ Get Region """
        if not self._region:
            region = tools.geometry.getRegion(self.area)
            self._region = region

        return self._region

    @property
    def date_range(self):
        """ Range of the date measure

        :rtype: ee.List
        """
        if not self._date_range:
            ordered = self.timeSeries.sort('system:time_start', True)
            def get_relative_date(img, ini):
                ini = ee.List(ini)
                date = img.date()
                measure = date.get(self.date_measure)
                return ini.add(measure)

            result = ordered.iterate(get_relative_date, ee.List([]))
            self._date_range = ee.List(result)

        return self._date_range


    @property
    def date_range_bitreader(self):
        """ Make a BitReader to encode the breakpoints

        Example:

        BitReader({0: 1999, 1: 2000, ... 12:2010}
        """
        if not self._date_range_bitreader:
            time_list = self.date_range.getInfo()
            reader_dict = {}
            for i, t in enumerate(time_list):
                key = '{}'.format(i)
                reader_dict[key] = {1: t}

            self._date_range_bitreader = bitreader.BitReader(reader_dict)

        return self._date_range_bitreader

    @property
    def core(self):
        """ Apply original LandTrendr Algorithm as given in GEE

        :return: Same as the original algorithm (An array image with one band
                 called 'LandTrendr' and the rest of the bands named:
                 bandname_fit (ej: B2_fit))
        """
        if not self._core:
            # Select the index band only (example: 'nbr') in the whole
            # collection
            timeserie_index = self.timeSeries.select(self.fit_band)

            img = ee.Algorithms.TemporalSegmentation.LandTrendr(
                timeSeries=timeserie_index,
                # self.timeSeries,
                maxSegments=self.maxSegments,
                spikeThreshold=self.spikeThreshold,
                vertexCountOvershoot=self.vertexCountOvershoot,
                preventOneYearRecovery=self.preventOneYearRecovery,
                recoveryThreshold=self.recoveryThreshold,
                pvalThreshold=self.pvalThreshold,
                bestModelProportion=self.bestModelProportion,
                minObservationsNeeded=self.minObservationsNeeded)

            self._core = img

        return self._core

    @property
    def breakdown(self):
        ''' This method breaks down the resulting array and returns a list of
        images

        returns an ImageCollection in which each image has the following bands:
            - year = image's year
            - {ind} = index (ndvi, evi, etc)
            - {ind}_fit = index ajustado
            - bkp = breakpoint (1: yes, 0: no).

        It assumes that each image in the collection represents only one year
        of the time series
        '''
        if not self._breakdown:
            core = self.core
            ltr = core.select('LandTrendr')
            rmse = core.select('rmse')
            n = self.timeSeries.size()
            ordered_ts = self.timeSeries.sort('system:time_start')
            ordered_list = ordered_ts.toList(n)
            seq = ee.List.sequence(0, n.subtract(1))
            def create(position, ini):
                ini = ee.List(ini)
                nextt = ee.Number(position).add(1)
                start = ee.Image.constant(ee.Number(position))
                end = ee.Image.constant(nextt)
                sli = ltr.arraySlice(1, start.toInt(), end.toInt(), 1)

                # CONVERT ARRAY TO IMG
                imgIndx = sli.arrayProject([0]).arrayFlatten(
                    [["year", self.fit_band, self.fit_band + "_fit", "bkp"]])

                date = ee.Image(ordered_list.get(position)).date().millis()
                result = imgIndx.addBands(rmse).set('system:time_start', date)
                return ini.add(result)
            collist = ee.List(seq.iterate(create, ee.List([])))
            col = ee.ImageCollection(collist)

            self._breakdown = col

        return self._breakdown

    @property
    def statistics(self):
        """ Compute statistics for this object """
        if not self._statistics:
            collection = self.slope
            self._statistics = statistics(collection, self.fit_band)

        return self._statistics

    @property
    def slope(self):
        """ Calculate slope of each segment in the LandTrendR fit

        Use:
        LandTrendr(imgcol, maxSegment, index, **kwargs).slope()

        Each image in the collection contains the following bands:
            - [0] index: original index
            - [1] {index}_fit: fitted index (example, ndvi_fit)
            - [2] bkp: breakpoint (1 change, 0 no change)
            - [3] year: image's year
            - [4] slope_before: slope of the segment 'before'. If *a1*
              is the main year, *slope_before = a1-a0*
            - [5] slope_after: slope of the segment 'after'. If *a1* is the
              main year, *slope_after = a2-a1*
            - [6] change: change's magnitud
              (between -1 (greater loose) and 1 (greater gain))
            - [7] angle: change's angle in radians
            - [8] rmse: root mean square error (original)

        :rtype: ee.ImageCollection

        """
        if not self._slope:
            breakdown = self.breakdown

            def add_date(img):
                """ Pass system:time_start property """
                date = img.get("system:time_start")
                return img.set("system:time_start", date)

            collection = breakdown.map(add_date)

            def compute(central, before=None, after=None):
                """ Compute slope with image before and after

                :param central: imagen 'central' (1)
                :type central: ee.Image

                :param before: imagen anterior (0)
                :type before: ee.Image

                :param after: imagen posterior (2)
                :type after: ee.Image

                :return: the central image with computed values:

                    **{band}**: original band value

                    **{band}_fit**: fitted band value

                    **slope_before**: slope to the image before

                    **slope_after**: slope to the image after

                    **change**: change magnitude (between -1 (greatest loss)
                    and 1 (greatest gain))

                    **angle**: change angle in radians
                :rtype: ee.Image
                """
                name_before = 'slope_before'
                name_after = 'slope_after'

                ### CENTRAL ################################################
                central_img = ee.Image(central)
                central_fit = central_img.select(self.fit_band + "_fit")
                central_year = central_img.select("year")
                central_date = central_img.get("system:time_start")

                #### BEFORE #################################################
                if before:
                    before_img = ee.Image(before)

                    before_fit = before_img.select(self.fit_band + "_fit")
                    before_year = before_img.select("year")

                    # difference (year) between central_img and before_img
                    before_diff_year = central_year.subtract(before_year)

                    # difference (fit) between central_img and before_img
                    before_diff_fit = central_fit.subtract(before_fit)\
                                                 .divide(before_diff_year)\
                                                 .select([0], [name_before])

                #### AFTER ##############################################
                if after:
                    after_img = ee.Image(after)

                    after_fit = after_img.select(self.fit_band + "_fit")
                    after_year = after_img.select("year")
                    after_diff_year = after_year.subtract(central_year)
                    after_diff_fit = after_fit.subtract(central_fit)\
                                              .divide(after_diff_year)\
                                              .select([0], [name_after])

                #### FIRST IMAGE #################
                if after and not before:
                    before_diff_fit = after_diff_fit.select([0], [name_before])

                elif before and not after:
                    after_diff_fit = before_diff_fit.select([0], [name_after])

                #### change ################################################
                pi = ee.Image.constant(math.pi)

                # angleS
                diff_before_radians = before_diff_fit.atan()
                diff_after_radians = after_diff_fit.atan()

                # Difference between after and before
                diff_fit = after_diff_fit.subtract(before_diff_fit)

                # Conditional mask
                gain_mask = diff_fit.gte(0)

                # each case images
                loss_img = diff_before_radians.subtract(diff_after_radians)\
                                              .subtract(pi)
                gain_img = pi.add(diff_before_radians)\
                             .subtract(diff_after_radians)

                # angle
                angle = ee.Image(loss_img.where(gain_mask, gain_img))\
                          .select([0], ["angle"])\
                          .toFloat()

                # Magnitude (between -1 y 1)
                gain_magnitude = pi.subtract(angle).divide(pi)
                loss_magnitude = pi.add(angle).divide(pi).multiply(-2)

                change = ee.Image(loss_magnitude.where(gain_mask, gain_magnitude))\
                           .select([0], ["change"])\
                           .toFloat()

                # gain
                pos_pi = ee.Number(3.14)
                gain = angle.lt(pos_pi).And(angle.gt(0)).toInt()

                # loss
                neg_pi = ee.Number(-3.14)
                loss = angle.gt(neg_pi).And(angle.lt(0)).toInt()

                img = central_img.addBands(before_diff_fit)\
                                 .addBands(after_diff_fit)\
                                 .addBands(change)\
                                 .addBands(angle)\
                                 .addBands(gain.select([0], ['gain']))\
                                 .addBands(loss.select([0], ['loss']))

                img = img.set("system:time_start", central_date)

                return img

            # Collection to List
            col_list = collection.toList(collection.size())

            new_col_list = ee.List([])

            # collection size
            tam = col_list.size().getInfo()  # - 1

            # first image
            im0 = col_list.get(0)

            # second image
            im1 = col_list.get(1)

            # compute first and second images
            im0 = compute(central=im0, after=im1)

            # add to list
            new_col_list = new_col_list.add(im0)

            # middle
            for im in range(1, tam - 1):
                im0 = col_list.get(im - 1)
                im1 = col_list.get(im)
                im2 = col_list.get(im + 1)

                im1 = compute(central=im1, before=im0, after=im2)
                new_col_list = new_col_list.add(im1)

            # last image
            im0 = col_list.get(tam - 2)
            im1 = col_list.get(tam - 1)

            im1 = compute(central=im1, before=im0)
            new_col_list = new_col_list.add(im1)

            newCol = ee.ImageCollection(new_col_list)

            self._slope = newCol

        return self._slope

    def total_bkp(self, collection=None):
        """ Compute the total number of breakpoint in each pixel

        :param collection: the collection that hold the breakpoint data, if
            None, it'll be computed.
        :type collection: ee.ImageCollection
        :return: A single Image with the total number of breakpoints in a band
            called `total_bkp`
        :rtype: ee.Image
        """
        col = collection if collection else self.slope
        sum_bkp = ee.Image(col.select("bkp").sum()).select([0], ["total_bkp"])
        return sum_bkp

    def stack(self):
        """ Generate an image holding the fitted index value for each year

        :return: An image with the following bands

            :{index}_{year}: fitted index value for that year
        :rtype: ee.Image
        """
        col = self.breakdown
        anio0 = ee.Date(col[0].get("system:time_start")).get("year").format()

        imgF = col[0].select([self.fit_band + "_fit"],
                             [ee.String(self.fit_band).cat(ee.String('_')).cat(anio0)])

        for i in range(1, len(col)):
            anio = ee.Date(col[i].get("system:time_start")).get(
                "year").format()
            img = col[i].select([self.fit_band + "_fit"],
                                [ee.String(self.fit_band).cat(ee.String('_')).cat(anio)])
            imgF = imgF.addBands(img)

        return imgF

    def breakpoints(self):
        """ Number of breakpoints

        Return an object with 2 properties:

        - image (ee.Image): ee.Image containing the following bands:

          - n_bkp: number of breakpoints in each pixel

        - total (ee.Number): maximum amount of breakpoints

        :rtype: namedtuple
        """
        # APLICO LANDTRENDR
        serie = self.breakdown
        col = ee.ImageCollection(serie)

        imgacum = ee.Image.constant(0).select([0], ["n_bkp"]).toUint8()
        # resta = ee.Image.constant(2)

        def check_bkp(img, acum):
            bkp = img.select("bkp")
            acum = ee.Image(acum)
            acum = acum.add(bkp)
            return acum

        n_bkp = ee.Image(col.iterate(check_bkp, imgacum))
        img_final = n_bkp#.subtract(resta)
        # geometry = ee.Image(self.timeSeries.first()).geometry().buffer(-5000)

        # workaround = ee.FeatureCollection(ee.List([ee.Feature(geometry)]))

        total = img_final.reduceRegion(ee.Reducer.max(),
                                       # geometry=geometry,
                                       geometry=self.area,
                                       scale=30,
                                       maxPixels=1e13).get("n_bkp")

        total = ee.Number(total).toInt8()

        bkps = namedtuple("Breakpoints", ["image", "total"])

        return bkps(image=img_final, total=total)

    def loss(self, skip_middle=True, skip_slope=0.05):
        """ Compute loss magnitude from a loss breakpoint until next
        breakpoint

        :param skip_middle: If the next date is also a loss, skip it.
        :type skip_middle: bool
        :param skip_slope: skip middle if the slope of the stretch before is
            greater than this value
        :type skip_slope: float
        :rtype: ee.ImageCollection
        """
        slope = self.slope.sort('system:time_start', True)
        last_date = ee.Number(self.date_range.get(-1)).toInt()
        fit_band = self.fit_band + '_fit'

        def over_slope(img, inidict):
            inidict = ee.Dictionary(inidict)
            ini = ee.List(inidict.get('images'))

            fitted = img.select([fit_band])
            year = img.select(['year'])

            # get image date to pass to the computed image
            time_start = img.date().millis()

            # statistics
            rmse = img.select(['rmse'])

            # date list from next to last
            date = img.date().get(self.date_measure).toInt()
            rrange = ee.List.sequence(date.add(1), last_date)

            # loss?
            loss_mask = img.select(['loss'])

            # loss before?
            loss_before = ee.Image(inidict.get('loss_before'))

            # update loss_before
            new_loss_before = loss_mask.select([0], ['loss_before'])

            times = tools.image.empty(0, ['times'])
            rest = tools.image.empty(0, ['next_bkp', 'loss'])

            initial = times.addBands(rest)

            # function to get the amount of loss and the elapsed loss time
            def over_range(d, ini):
                # cast
                ini = ee.Image(ini)
                d = ee.Number(d)

                date_range = ee.Image(d.subtract(date)).select([0], ['next_bkp'])

                # get slope image for this d (date)
                img_s = ee.Image(
                    slope.filter(
                        ee.Filter.calendarRange(d, d, self.date_measure)
                    ).first())

                # get image date before
                img_before = ee.Image(
                    slope.filter(
                        ee.Filter.calendarRange(d.subtract(1),
                                                d.subtract(1),
                                                self.date_measure)
                    ).first())

                img_before_loss = img_before.select(['loss'])
                img_before_bkp = img_before.select(['bkp'])

                # retain only loss pixels
                img_s = img_s.updateMask(loss_mask)

                # retain pixels that are not preceeded by loss
                # img_slope_before = img_s.select('slope_before')
                mask = loss_before.Not()#.And(img_slope_before.gte(skip_slope))
                img_s = img_s.updateMask(mask)

                # catch loss
                this_loss_mask = img_s.select(['loss'])

                # catch gain
                this_gain_mask = img_s.select(['gain'])

                # catch breakpoint
                bkp = img_s.select(['bkp'])

                # make a mask only if it is the first occurrence of gain
                first_time = ini.select(['times']).eq(0).toInt()

                # CASE 1
                case1 = this_gain_mask.And(img_before_loss)

                # CASE 2
                case2 = this_loss_mask.And(img_before_bkp.Not())

                # CASE 3
                case3 = this_loss_mask.Not().And(this_gain_mask.Not()).And(bkp)

                # EXTRA CASE
                if not skip_middle:
                    extra_case = this_loss_mask.And(img_before_loss)
                    retain = case1.Or(case2).Or(case3).Or(extra_case).And(first_time)
                else:
                    # retain
                    retain = case1.Or(case2).Or(case3).And(first_time)

                # update times
                times = ini.select(['times']).add(retain)

                # get loss magnitude
                new_fitted = img_s.select([fit_band])
                loss = fitted.subtract(new_fitted) \
                    .select([0], ['loss_magnitude']) \
                    .updateMask(retain).unmask()

                range_i = date_range.updateMask(retain).unmask().toInt()

                to_add = times.addBands(range_i).addBands(loss)

                final = ini.select(['times', 'next_bkp', 'loss']).add(to_add)

                return final

            next_loss = ee.Image(rrange.iterate(over_range, initial)) \
                          .addBands(year)\
                          .addBands(rmse)\
                          .addBands(loss_before)\
                          .set('system:time_start', time_start)

            # set images
            inidict = inidict.set('images', ini.add(next_loss))

            return ee.Dictionary(inidict).set('loss_before', new_loss_before)

        inidict = ee.Dictionary({
            'images': ee.List([]),
            'loss_before': tools.image.empty(0, ['loss_before'])
        })

        newdict = ee.Dictionary(slope.iterate(over_slope, inidict))
        images = ee.List(newdict.get('images'))
        return ee.ImageCollection.fromImages(images)

    # ENCODED IMAGES

    def breakpoints_image(self):
        ''' Create an image with breakpoints occurrence encoded by a BitReader
        '''
        breakdown = self.breakdown
        encoder = self.date_range_bitreader
        time_list = self.date_range.getInfo()

        # trim last and first value of time_list because are bkps always
        time_list = time_list[1:-1]

        # Create a dict with encoded values for each time
        values = {}
        for t in time_list:
            encoded = encoder.encode(t)
            values[t] = encoded

        valuesEE = ee.Dictionary(values)

        ini_img = tools.image.empty(from_dict={'bkp':0})

        def over_col(img, ini):
            ini = ee.Image(ini)
            date = img.date().get(self.date_measure).format()
            bkp = img.select(['bkp'])
            date_value = valuesEE.get(date)
            date_value_img = ee.Image.constant(date_value).select([0], ['bkp'])
            masked_date_value = date_value_img.updateMask(bkp).unmask()

            return ini.add(masked_date_value)

        # slice first and last of breakdown because bkps are always 1
        breakdown = ee.ImageCollection.fromImages(
                       breakdown.toList(breakdown.size())\
                             .slice(1, -1))

        result = ee.Image(breakdown.iterate(over_col, ini_img))
        result = result.clip(ee.Geometry.Polygon(self.region))

        return result

    def lossdate_image(self, threshold, bandname='lossyear', skip_middle=True):
        """ make an encoded loss image holding dates of disturbance """
        # slope = self.slope
        loss = self.loss(skip_middle)

        encoder = self.date_range_bitreader
        time_list = self.date_range.getInfo()

        # Create a dict with encoded values for each time
        values = {}
        for t in time_list:
            encoded = encoder.encode(t)
            values[t] = encoded

        valuesEE = ee.Dictionary(values)

        def over_col(img, ini):
            # last lossyear img
            ini = ee.Image(ini)

            # get date (example: 2010)
            date = img.date().get(self.date_measure).format()
            date_value = valuesEE.get(date)
            date_value_img = ee.Image.constant(date_value) \
                .select([0], [bandname])

            condition = img.select(['loss']).gt(threshold)

            masked_date_value = date_value_img.updateMask(condition).unmask()

            return ini.add(masked_date_value)

        initial = tools.image.empty(0, [bandname])
        result = ee.Image(loss.iterate(over_col, initial))
        result = result.clip(ee.Geometry.Polygon(self.region))

        return result