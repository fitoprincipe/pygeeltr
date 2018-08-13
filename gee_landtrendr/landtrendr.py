#!/usr/bin/env python
# -*- coding: utf-8 -*-

from geetools import tools
import math
from collections import namedtuple
from itertools import chain

import ee

# Initialize EE if not initialized
if not ee.data._initialized: ee.Initialize()


def statistics(collection, band=None, suffix='_fit'):
    """ Compute sum of squares and R2 for values of every pixel in a
    ImageCollection that has already be fitted

    :param collection: collection
    :type collection: ee.ImageCollection
    :param band: name of the band that has been fitted. If None, the first
        band will be used
    :type band: str
    :param suffix: suffix of band that holds the fitted values
    :type suffix: str
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

    # Mean
    mean = justbands.mean()

    # SSTOT
    sstot_image = tools.empty_image(0, ['sstot'])

    def sstot(img, sstoti):
        sstoti = ee.Image(sstoti) # cast
        mean_i = mean.select([band])
        val = img.select([band])
        compute = val.subtract(mean_i).pow(2).select([0], ['sstot'])
        return sstoti.add(compute)

    sstot_image = ee.Image(justbands.iterate(sstot, sstot_image))

    # SSREG
    ssreg_image = tools.empty_image(0, ['ssreg'])

    def ssreg(img, ssregi):
        ssregi = ee.Image(ssregi) # cast
        mean_i = mean.select([band])
        val = img.select([fitted_band])
        compute = val.subtract(mean_i).pow(2).select([0], ['ssreg'])
        return ssregi.add(compute)

    ssreg_image = ee.Image(justbands.iterate(ssreg, ssreg_image))

    # SSRES
    ssres_image = tools.empty_image(0, ['ssres'])

    def ssres(img, ssresi):
        ssresi = ee.Image(ssresi) # cast
        val = img.select([band])
        fit = img.select([fitted_band])
        compute = val.subtract(fit).pow(2).select([0], ['ssres'])
        return ssresi.add(compute)

    ssres_image = ee.Image(justbands.iterate(ssres, ssres_image))

    # DIVIDE
    division = ssres_image.divide(sstot_image).select([0], ['r2'])

    # 1-division
    r2 = tools.empty_image(1, ['r2']).subtract(division).toFloat()

    result = namedtuple("Statistics", ["ssres", "r2"])

    return result(ssres_image, r2)


class LandTrendr(object):
    """
    LandTrendr (Landsat-based detection of Trends in Disturbance and
    Recovery) is an algorithm that temporally segments a time-series of
    images by extracting the spectral trajectories of change over time.
    The first band of each image is used to find breakpoints, and those
    breakpoints are used to perform fitting on all subsequent bands.
    The breakpoints are returned as a 2-D matrix of 4 rows and as many
    columns as images. The first two rows are the orignial X and Y
    values. The third row contains the Y values fitted to the estimated
    segments, and the 4th row contains a 1 if the corresponding point
    was used as a segment vertex or 0 if not.

      **Parameters:**

    - timeSeries (ImageCollection):

      Collection from which to extract trends (it's assumed thateach
      image in the collection represents one year). The first band is
      usedto find breakpoints, and all subsequent bands are fitted
      using those breakpoints.

    - index (str):

      band to use for regression

      **Optionals:**

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

    def __init__(self, timeserie, index, area=None, **kwargs):

        self.timeSeries = timeserie
        self.index = index

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
        # self.col_size = tools.execli(self.timeSeries.size().getInfo)()

        if area:
            self.area = area
        else:
            self.area = ee.Image(self.timeSeries.first()).geometry()

        # self.region = tools.execli(self.area.getInfo)()["coordinates"]
        self.region = tools.getRegion(self.area)

    @classmethod
    def Principe(cls, timeserie, index, area=None):
        """ Factory Method para crear una clase con los parametros
        definidos por Principe
        """
        newobj = cls(timeserie, index, area,
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
    def Kennedy(cls, timeserie, index, area=None):
        """ Factory Method para crear una clase con los parametros
        definidos por Kennedy (original)
        """
        newobj = cls(timeserie, index, area,
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
    def Liang(cls, timeserie, index, area=None):
        """ Factory Method para crear una clase con los parametros
        definidos por Liang
        """
        newobj = cls(timeserie, index, area,
                     maxSegments=4,
                     spikeThreshold=0.9,
                     vertexCountOvershoot=0,
                     preventOneYearRecovery=False,
                     recoveryThreshold=0.25,
                     pvalThreshold=0.1,
                     bestModelProportion=0.75,
                     minObservationsNeeded=6)

        return newobj

    def CORE(self):
        """ Apply original LandTrendr Algorithm as given in GEE
        :return: namedtuple

        :result: Same as the original algorithm (An array image with one band
                 called 'LandTrendr' and the rest of the bands named:
                 bandname_fit (ej: B2_fit))

        :year_list: list of millisecond dates (ee.List)
        """
        # if self.col_size > 1:  # 1 Request

        # Select the index band only (example: 'nbr') in the whole collection
        timeserie_index = self.timeSeries.select(self.index)

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

        '''
        def anioslist(img, listanios):
            time = ee.Date(img.get("system:time_start")).millis()
            return ee.List(listanios).add(time)

        anios = ee.List(self.timeSeries.iterate(anioslist, ee.List([])))

        core = namedtuple("CORE", ["result", "year_list"])

        return core(result=img, year_list=anios)
        # else:
        #     raise ValueError("The time serie must have more than 1 image")
        '''
        return img

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
        # core = self.CORE().result
        core = self.CORE()
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
                [["year", self.index, self.index + "_fit", "bkp"]])

            date = ee.Image(ordered_list.get(position)).date().millis()
            result = imgIndx.addBands(rmse).set('system:time_start', date)
            return ini.add(result)
        collist = ee.List(seq.iterate(create, ee.List([])))
        col = ee.ImageCollection(collist)

        return col

    def statistics(self):
        """ Compute statistics for this object """
        collection = self.slope()
        return statistics(collection, self.index)

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

        :rtype: ee.ImageCollection

        """

        # APLICO LANDTRENDR
        #listaImgsPy = self.breakdown()
        #listaImgs = ee.List(listaImgsPy)  # ee.List
        #colLTR = ee.ImageCollection(listaImgs)  # ee.ImageCollection
        colLTR = self.breakdown()

        def addTime(img):
            """ Pass system:time_start property """
            prop_anio = img.get("system:time_start")
            # factor = ee.Image.constant(1000)
            # img = img.select(self.index, self.index + "_fit", "bkp") #.addBands(anioI)  # .multiply(factor).toInt()

            return img.set("system:time_start", prop_anio)

        colInt = colLTR.map(addTime)

        # DEBUG
        # ultima = ee.Image(colInt.toList(50).get(colInt.size().getInfo()-1))
        # print "fecha de la ultima img:", ee.Date(ultima.date()).get(
        # "year").getInfo()

        # NOMBRES PARA LAS BANDAS
        adelante = ee.String("slope_after")
        atras = ee.String("slope_before")

        def calcIndices(**kwargs):
            """ Función para calcular los indices de los slope comparando
            los indices de la imagen anterior y posterior a la central. Si no
            hay anterior o posterior ese segmento lo pone en 0

            :param central: imagen 'central' (1)
            :type central: ee.Image

            :param ant: imagen anterior (0)
            :type ant: ee.Image

            :param pos: imagen posterior (2)
            :type pos: ee.Image

            :return: la imagen central con los indices calculados:

                **index**: index original

                **indice_fit**: index ajustado

                **slope_before**: diferencia atras

                **slope_after**: diferencia adelante

                **change**: magnitud del change (entre -1 (> perdida)
                y 1 (> ganancia))

                **angle**: angle del change en radianes
            :rtype: ee.Image
            """

            ### CENTRAL ################################################
            img1orig = ee.Image(kwargs.get("central"))
            img1fit = img1orig.select(self.index + "_fit")
            anio1 = img1orig.select("year")
            d1 = img1orig.get("system:time_start")

            #### ATRAS #################################################
            # if kwargs.has_key("ant"):
            if "ant" in kwargs.keys():
                img0 = ee.Image(kwargs.get("ant"))

                img0fit = img0.select(self.index + "_fit")
                anio0 = img0.select("year")

                # CALCULO LA DIFERENCIA DE AÑOS ENTRE IM1 E IM0
                dif_anio_At = anio1.subtract(anio0)

                # CALCULO LA DIFERENCIA DEL INDICE AJUSTADO
                difAt = (img1fit.subtract(img0fit)
                         .divide(dif_anio_At).select([0], [atras]))

            #### ADELANTE ##############################################
            # if kwargs.has_key("pos"):
            if "pos" in kwargs.keys():
                img2 = ee.Image(kwargs.get("pos"))

                img2fit = img2.select(self.index + "_fit")
                anio2 = img2.select("year")
                dif_anio_Ad = anio2.subtract(anio1)
                difAd = (img2fit.subtract(img1fit)
                         .divide(dif_anio_Ad).select([0], [adelante]))

            #### SI NO HAY IMG ANTERIOR Y SI POSTERIOR #################
            # if not kwargs.has_key("ant") and kwargs.has_key("pos"):
            if "ant" not in kwargs.keys() and "pos" in kwargs.keys():
                difAt = difAd.select([0], [atras])
            # elif not kwargs.has_key("pos") and kwargs.has_key("ant"):
            elif "pos" not in kwargs.keys() and "ant" in kwargs.keys():
                difAd = difAt.select([0], [adelante])

            #### change ################################################
            pi = ee.Image.constant(math.pi)

            # angleS
            difAtRad = difAt.atan()
            difAdRad = difAd.atan()

            # DIFERENCIA ATRAS Y ADELANTE
            dif = difAd.subtract(difAt)

            # MASCARAS CONDICIONALES
            cond_ganancia = dif.gte(0)

            # IMAGENES PARA CADA CASO
            img_perdida = difAtRad.subtract(difAdRad).subtract(pi)
            img_ganancia = pi.add(difAtRad).subtract(difAdRad)

            # angle
            angle = (ee.Image(img_perdida.where(cond_ganancia, img_ganancia))
                      .select([0], ["angle"])
                      .toFloat())

            # MAGNITUD (entre -1 y 1)
            mag_ganancia = pi.subtract(angle).divide(pi)
            mag_perdida = pi.add(angle).divide(pi).multiply(-2)

            change = (ee.Image(mag_perdida.where(cond_ganancia, mag_ganancia))
                      .select([0], ["change"])
                      .toFloat())

            img = img1orig.addBands(difAt).addBands(difAd).addBands(
                change).addBands(angle)
            img = img.set("system:time_start", d1)

            return img

        # CONVIERTO LA COL A UNA LISTA
        colListE = colInt.toList(50)

        # CREO UNA LISTA VACIA PARA IR AGREGANDO LAS IMG RESULTADO Y CREAR
        # LA NUEVA COLECCION
        newColList = ee.List([])

        # TAMAÑO DE LA COLECCION
        tam = tools.execli(colListE.size().getInfo)()  # - 1

        # PRIMER IMAGEN
        im0 = colListE.get(0)

        # SEGUNDA IMAGEM
        im1 = colListE.get(1)

        # CALCULA LOS INDICES ENTRE LA PRIMERA Y LA SEGUNDA IMG
        im0 = calcIndices(central=im0, pos=im1)

        # AGREGA EL RESULTADO A LA LISTA VACIA
        newColList = newColList.add(im0)

        # MEDIO
        for im in range(1, tam - 1):
            im0 = colListE.get(im - 1)
            im1 = colListE.get(im)
            im2 = colListE.get(im + 1)

            im1 = calcIndices(central=im1, ant=im0, pos=im2)
            newColList = newColList.add(im1)

        # ULTIMA IMAGEN
        im0 = colListE.get(tam - 2)
        im1 = colListE.get(tam - 1)

        im1 = calcIndices(central=im1, ant=im0)
        newColList = newColList.add(im1)

        newCol = ee.ImageCollection(newColList)

        return newCol

    def total_bkp(self, collection=None):
        """ Compute the total number of breakpoint in each pixel

        :param collection: the collection that hold the breakpoint data, if
            None, it'll be computed.
        :param collection: ee.ImageCollection
        :return: A single Image with the total number of breakpoints in a band
            called `total_bkp`
        :rtype: ee.Image
        """
        col = collection if collection else self.slope()
        sum_bkp = ee.Image(col.select("bkp").sum()).select([0], ["total_bkp"])
        return sum_bkp

    def stack(self):
        """ Generate an image holding the fitted index value for each year

        :return: An image with the following bands

            :{index}_{year}: fitted index value for that year
        """
        col = self.breakdown()
        anio0 = ee.Date(col[0].get("system:time_start")).get("year").format()

        imgF = col[0].select([self.index + "_fit"],
                             [ee.String(self.index).cat(ee.String('_')).cat(anio0)])

        for i in range(1, len(col)):
            anio = ee.Date(col[i].get("system:time_start")).get(
                "year").format()
            img = col[i].select([self.index + "_fit"],
                                [ee.String(self.index).cat(ee.String('_')).cat(anio)])
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
        serie = self.breakdown()
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

        # total = tools.execli(ee.Number(total).getInfo)()#.add(ee.Number(2)).getInfo()
        # total = int(total)

        bkps = namedtuple("Breakpoints", ["image", "total"])

        return bkps(image=img_final, total=total)

    def break2band(self):
        """ Generate an image holding the data for each found breakpoint

        Return an Image with the following bands:

        - year_{breakpoint}: year for each breakpoint. Number of breakpoints in
          the whole image may differ from the number of breakpoints in
          each pixel. Therefore, bands will correspond to breakpoint in
          the whole image.

            Example:

            - bkps in the whole image: 5
            - bkps in one pixel: 3
            - values of bands in that pixel:

              - year_1: 0
              - year_2: 0
              - year_3: 1999
              - year_4: 2005
              - year_5: 2017

        - change_{id}: change value for each breakpoint
        - backward_{id}: backward change value
        - forward_{id}: forward change value

        :rtype: ee.Image
        """
        # APLICO LANDTRENDR
        # breakdown = self.breakdown()
        # col = ee.ImageCollection(serie)

        sufijos = ["year", "change", "backward", "forward", "fit"]
        sufijos = ee.List(sufijos)

        col = self.slope()

        # CONVIERTO LA COL A UNA LISTA
        colL = col.toList(100)

        # PRIMER IMG DE LA LISTA
        img0 = ee.Image(colL.get(0))

        # IMAGEN NUEVA CON LA BANDA EN LA QUE SE ACUMULARAN LOS bkp
        img_ini = ee.Image(1).select([0],["id_bkp"]).toUint8()

        # LISTA CON UNA IMG: LA PRIMER IMAGEN ORIGINAL SUMADA LA BANDA id_bkp
        acum = ee.List([img0.addBands(img_ini)])

        # FUNCION PARA ACUMULAR bkp
        def addId(img, acum):
            # CAST acum
            listacum = ee.List(acum)

            # OBTENGO LA ULTIMA IMG DE LA LISTA
            ultima = ee.Image(listacum.get(-1))

            # OBTENGO LA BANDA id_bkp DE LA ULTIMA IMG
            bkp_ultima = ultima.select("id_bkp")

            # OBTENGO LA BANDA bkp DE LA IMG ACTUAL (0 Y 1) Y LA RENOMBRO
            bkp_img = ee.Image(img).select(["bkp"], ["id_bkp"]).toUint8()

            # SUMO LA BANDA DE BREAKPOINTS ACUMULADOS CON LA DE BREAKPOINT ACTUAL
            nueva_bkp = ee.Image(bkp_img.add(bkp_ultima)).toUint8()

            # AGREGO LA BANDA SUMADA A LA IMG ACTUAL
            nueva_img = ee.Image(img).addBands(nueva_bkp)

            # DEVUELVO LA LISTA CON LA IMAGEN AGREGADA
            return listacum.add(nueva_img)

        colL = ee.List(colL.slice(1).iterate(addId, acum))

        col = ee.ImageCollection(colL)

        # OBTENGO LA CANT MAXIMA DE BKPS
        maxbkp = self.breakpoints().total  # ee.Number

        # LISTA DE INDICES
        # indices = range(0, maxbkp)
        # indices_str = [str(i) for i in indices]
        def int2str(integer):
            return ee.Number(integer).toInt().format()
        indices_str = ee.List.sequence(0, maxbkp.subtract(1)).map(int2str)

        # CONVIERTO EN 0 A LOS VALORES EN LOS PIXELES QUE bkp = 0
        def conversion(img):
            condicion = img.select("bkp").eq(0)
            return img.where(condicion, ee.Image.constant(0))

        col = col.map(conversion)

        # CONVIERTO LA COL A UN ARRAY
        array = col.toArray()

        ejeImgs = 0
        ejeBandas = 1

        # ORDENO EL ARRAY SEGUN LA BANDA id_bkp
        bkps = array.arraySlice(ejeBandas, 2, 3)
        ordenado_bkp = array.arraySort(bkps)

        # IMG DE LONGITUD TOTAL DEL ARRAY (CANT DE IMGS)
        longitud = ordenado_bkp.arrayLength(ejeImgs)

        # CORTO EL ARRAY DEJANDO SOLO LOS AÑOS INDICADOS en maxbkp
        cortado_bkp = ordenado_bkp.arraySlice(ejeImgs,
                                              longitud.subtract(maxbkp),
                                              longitud)

        # DEBUG
        # print(cortado_bkp.getInfo())

        # ORDENO EL ARRAY CORTADO SEGUN LOS AÑOS
        anios = cortado_bkp.arraySlice(ejeBandas, 3, 4)
        ordenado_anios = cortado_bkp.arraySort(anios).arraySlice(ejeBandas, 3, 4)

        # CORTO LA BANDA DE CAMBIO
        cambio = cortado_bkp.arraySlice(ejeBandas, 6, 7)
        array_final = ordenado_anios.arrayCat(cambio, ejeBandas)

        # CORTO LAS BANDAS slope_before Y slope_after
        slope_before = cortado_bkp.arraySlice(ejeBandas, 4, 5)
        slope_after = cortado_bkp.arraySlice(ejeBandas, 5, 6)
        fit = cortado_bkp.arraySlice(ejeBandas, 1, 2)
        array_final = array_final.arrayCat(slope_before, ejeBandas)
        array_final = array_final.arrayCat(slope_after, ejeBandas)
        array_final = array_final.arrayCat(fit, ejeBandas)

        # TRANSFORMO EN UNA IMAGES
        # img = ordenado_anios.arrayFlatten([indices_str,["a"]])
        # img = ordenado_anios.arrayTranspose().arrayFlatten([["a"], indices_str])
        img = array_final.arrayTranspose().arrayFlatten([sufijos, indices_str])
        return img

    def stretches(self, min_threshold=0.05, max_threshold=0.2):
        """ This method characterize the segmentation stretches and returns as
        many images as the amount of stretches. Categories are:

        * 1: no change: -min_threshold <= value <= min_threshold
        * 2: soft loss: -max_threshold <= value < -min_threshold
        * 3: steep loss: value < -max_threshold
        * 4: soft gain: min_threshold < value <= max_threshold
        * 5: steep gain: max_threshold < value

        :param min_threshold: divides 'no change' and 'soft change'
        :type min_threshold: float
        :param max_threshold: divides 'soft change' and 'steep change'
        :type max_threshold: float

        Return a new class called 'Stretches' with the following properties:

        - img_list: a list of images containing the stretches.

          Each image has the following bands:

            - t{n}_slope: slope of stretch n
            - t{n}_duration: duration of stretch n (in years)
            - t{n}_cat: category for stretch n

        - image: an image with unified results. Will have as many bands as
          found stretches, times 3. In case stretches in one pixel are
          less than stretches in the whole image, the last will be empty.

          Example:

          The whole image has 4 stretches. The pixel that has 2
          stretches will have the following values:

            - t1_slope: value
            - t2_slope: value
            - t3_slope: 0
            - t4_slope: 0

        :rtype: namedtuple
        """
        # Threshold images
        umb1pos = ee.Image.constant(min_threshold).toFloat()
        umb2pos = ee.Image.constant(max_threshold).toFloat()

        umb1neg = ee.Image.constant(-min_threshold).toFloat()
        umb2neg = ee.Image.constant(-max_threshold).toFloat()

        # Number of breakpoints in each pixel
        bkps = self.breakpoints()
        img_bkp = bkps.image

        # Total number of breakpoints
        total_bkp = bkps.total.getInfo()
        max_tramos = total_bkp - 1

        bandas_a = ["year_"+str(i) for i in range(total_bkp)]
        bandas_fit = ["fit_"+str(i) for i in range(total_bkp)]
        bandas = bandas_a + bandas_fit

        b2b = self.break2band().select(bandas)

        imagen = b2b.addBands(img_bkp)

        # OBTENGO TANTAS IMAGENES COMO TRAMOS HAYA EN LA COLECCION
        listaimg = []

        for t, bkp in enumerate(range(2, total_bkp+1)):
            tramos = t+1
            mask = imagen.select("n_bkp").eq(bkp)
            masked = imagen.updateMask(mask)
            img_tramos = ee.Image.constant(0)
            for id, tramo in enumerate(range(tramos, 0, -1)):
                # NOMBRE DE LA BANDA
                nombre_tramo = "t"+str(id+1)

                # INDICES INI Y FIN
                ini = max_tramos-tramo
                fin = ini + 1

                # SELECCIONO LAS BANDAS
                a_ini = ee.Image(masked.select("year_"+str(ini))).toFloat()
                a_fin = ee.Image(masked.select("year_"+str(fin))).toFloat()
                fit_ini = masked.select("fit_"+str(ini))
                fit_fin = masked.select("fit_"+str(fin))

                # COMPUTO
                lapso = a_fin.subtract(a_ini)
                dif_total = fit_fin.subtract(fit_ini)

                dif_anual = dif_total.divide(lapso)
                slope = dif_anual.select([0],[nombre_tramo+"_slope"])

                duracion = ee.Image(lapso.select([0],[nombre_tramo+"_duration"])).toUint8()

                # CATEGORIZACION
                uno = slope.gte(umb1neg).And(slope.lte(umb1pos))
                dos = slope.lt(umb1neg).And(slope.gte(umb2neg))
                tres = slope.lt(umb2neg)
                cuatro = slope.gt(umb1pos).And(slope.lte(umb2pos))
                cinco = slope.gt(umb2pos)

                cat = uno.add(
                      dos.multiply(2)).add(
                      tres.multiply(3)).add(
                      cuatro.multiply(4)).add(
                      cinco.multiply(5))

                cat = ee.Image(cat).select([0],[nombre_tramo+"_cat"]).toUint8()

                img_tramos = img_tramos.addBands(slope).addBands(duracion).addBands(cat)

            # SACO LA PRIMER BANDA QUE SE CREA ANTES DE INICIAR EL LOOP
            bandas = img_tramos.bandNames().slice(1)
            img_tramos = img_tramos.select(bandas)

            listaimg.append(img_tramos)

        # CREO UNA IMAGEN VACIA CON TODAS LAS BANDAS
        bandas_cat = ["t{0}_cat".format(str(n+1)) for n in range(max_tramos)]
        bandas_slope = ["t{0}_slope".format(str(n+1)) for n in range(max_tramos)]
        bandas_duracion = ["t{0}_duration".format(str(n+1)) for n in range(max_tramos)]

        bandas_todas = zip(bandas_slope, bandas_duracion, bandas_cat)

        bandas_todas = list(chain.from_iterable(bandas_todas))

        # bandas_todas = bandas_cat+bandas_slope+bandas_duracion

        lista_ini = ee.List(bandas_todas).slice(1)

        img_ini = ee.Image.constant(0).select([0],["t1_slope"])

        def addB(item, ini):
            ini = ee.Image(ini)
            img = ee.Image.constant(0).select([0],[item])
            return ini.addBands(img)

        img_final = ee.Image(lista_ini.iterate(addB, img_ini))

        for tramos, img in enumerate(listaimg):
            tramos += 1
            faltan = max_tramos - tramos

            # print "faltan", faltan
            # nombre = "falta{}tramo".format(str(faltan))

            if faltan == 0:
                img = tools.mask2zero(img)
                img_final = img_final.add(img)
                # funciones.asset(img, nombre, "users/rprincipe/Pruebas/"+nombre, self.region)
                continue

            for tramo in range(faltan):
                tramo += 1
                i = tramo + tramos
                bcat = ee.Image.constant(0).select([0],["t{0}_cat".format(str(i))])
                bslope = ee.Image.constant(0).select([0],["t{0}_slope".format(str(i))])
                bdur = ee.Image.constant(0).select([0],["t{0}_duration".format(str(i))])
                total = bcat.addBands(bslope).addBands(bdur)
                img = img.addBands(total)

            img = tools.mask2zero(img)

            # funciones.asset(img, nombre, "users/rprincipe/Pruebas/"+nombre, self.region)

            # bandas = img.getInfo()["bands"]
            # print "tramos:", tramos, [i["id"] for i in bandas]

            img_final = img_final.add(img)

        # return

        resultado = namedtuple("Stretches", ["img_list", "image"])

        return resultado(listaimg, img_final)
