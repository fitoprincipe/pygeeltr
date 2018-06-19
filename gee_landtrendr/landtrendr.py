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

    Argumentos:

    OBLIGATORIOS:

    timeSeries (ImageCollection):
    Collection from which to extract trends (it's assumed thateach
    image in the collection represents one year). The first band is
    usedto find breakpoints, and all subsequent bands are fitted
    using those breakpoints.

    index (str):
    nombre del index que se utiliza para hacer la regresion

    OPCIONALES:
    maxSegments (Integer):
    Maximum number of segments to be fitted on the time series.

    spikeThreshold (Float, default: 0.9):
    Threshold for dampening the spikes (1.0 means no dampening).

    vertexCountOvershoot (Integer, default: 3):
    The inital model can overshoot the maxSegments + 1 vertices by this
    amount. Later, it will be prunned down to maxSegments + 1.

    preventOneYearRecovery (Boolean, default: false):
    Prevent segments that represent one year recoveries.

    recoveryThreshold (Float, default: 0.25):
    If a segment has a recovery rate faster than 1/recoveryThreshold
    (in years), then the segment is disallowed.

    pvalThreshold (Float, default: 0.1):
    If the p-value of the fitted model exceeds this threshold, then the
    current model is discarded and another one is fitted using the
    Levenberg-Marquardt optimizer.

    bestModelProportion (Float, default: 1.25):
    Takes the model with most vertices that has a p-value that is at
    most this proportion away from the model with lowest p-value.

    minObservationsNeeded (Integer, default: 6):
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

        # PROPIEDADES CONSTANTES
        self.year_axis = 1
        self.band_axis = 0
        self.col_size = tools.execli(self.timeSeries.size().getInfo)()

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
        if self.col_size > 1:  # 1 Request

            # Select the index band only (example: 'nbr') in the whole collection
            timeserie_index = self.timeSeries.select(self.index)

            img = ee.Algorithms.Test.LandTrendr(
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

            listaTimeStart = ee.List([])

            def anioslist(img, listanios):
                time = ee.Date(img.get("system:time_start")).millis()
                return ee.List(listanios).add(time)

            anios = ee.List(self.timeSeries.iterate(anioslist, listaTimeStart))

            core = namedtuple("CORE", ["result", "year_list"])

            return core(result=img, year_list=anios)
        else:
            raise ValueError("The time serie must have more than 1 image")

    def breakdown(self):
        """ This method breaks down the resulting array and returns a list of
        images

        :return: A python list of images with the given bands:

            - year = image's year

            - {ind} = index (ndvi, evi, etc)

            - {ind}_fit = index ajustado

            - bkp = breakpoint (1: yes, 0: no).

            NO 5+. banda_aj = (nombre de la banda) ajustado. ej: SWIR_aj
        :rtype: list
        """

        # Get result and year list from the core method
        img = self.CORE().result
        years = self.CORE().year_list

        if img is not None:

            # band names
            bandas = img.bandNames()

            # OPCION 2
            # CREA UNA LISTA CON EL NOMBRE DE LAS BANDAS
            infobandas = tools.execli(bandas.getInfo)()
            # bandasN = funciones.decode_list(infobandas)  # 2 req
            bandasN = infobandas

            # SELECT LANDTRENDR BAND

            # trend = img.select([bandas.get(0)]) # ESTO FUNCIONABA
            trend = img.select(bandasN[0])  # ESTO FUNCIONA

            # PARA CADA IMG DE LA COLECCION

            listaImgs = []

            for indx in range(0, self.col_size):
                # CORTO EL ARRAY EN EL EJE anios, EL AÑO = indx
                imgIndx = trend.arraySlice(self.year_axis, indx, indx + 1)

                # CONVIERTO EL ARRAY A UNA IMG CON BANDAS
                imgIndx = imgIndx.arrayProject([0]).arrayFlatten(
                    [["year", self.index, self.index + "_fit", "bkp"]])

                # SUMA UNO A LA BANDA anio PORQUE LANDTRENDR RESTA UNO
                # IMPORTANTE! : LA BANDA anio PASA AL FINAL
                suma = ee.Image.constant(1).select([0],["year"])
                newanio = imgIndx.select("year").add(suma)

                imgIndx = tools.replace(imgIndx, "year", newanio)

                # PARA CADA BANDA DE LA IMG OBTENGO LOS VALORES AJUSTADOS

                sts = ee.Date(years.get(indx))

                # LE AGREGO LA PROPIEDAD system:time_start
                # print "238. system:time_start", date
                imgIndx = imgIndx.set("system:time_start", sts.millis())

                # print "miliseg", sts.millis().getInfo()
                # print "año", anios.get(indx)

                # UNA VEZ Q TENGO LA IMG FINAL LA AGREGO A UNA LISTA
                listaImgs.append(imgIndx)

            return listaImgs

    def extract_values_at(self, point):
        """ Dado un *punto* este metodo extre los valores del index, index
        ajustado y año de toda la serie temporal

        :param point: el punto del cual se quiere extraer la informacion
        :type point: ee.Geometry.Point

        :return: listaVal, listaValFit, listaAnios

            **listaVal:** lista de valores del index

            **listaValFit:** lista de valores del index ajustado por LandTrendr

            **listaAnios:** lista de años
        :rtype: tuple
        """
        # APLICO LANDTRENDR
        # img, anios = self.CORE()
        img = self.CORE().result
        anios = self.CORE().year_list

        # EXTRAIGO LA INFO DEL PUNTO
        extracto = tools.execli(img.reduceRegion(ee.Reducer.first(), point, 30).getInfo)()

        indice = extracto["LandTrendr"][1]
        indice_fit = extracto["LandTrendr"][2]
        anios = extracto["LandTrendr"][0]
        aniosPy = []

        # AUMENTO EL AÑO EN 1
        for anio in anios:
            aniosPy.append(anio + 1)

        return indice, indice_fit, aniosPy

    def point_to_dict(self, point):
        ''' '''
        pass

    def statistics(self):
        collection = self.slope()
        return statistics(collection, self.index)

    def statistics_array(self):
        """ Residual sum of square and coefficient of determination over
        LandTrendR fit curves

        :return:
            :ssres: residual sum of square in an ee.Image
            :r2: coefficient of determination (r2)
        :rtype: namedtuple
        """
        # APLICO LANDTRENDR
        # img, anios = self.CORE()
        img = self.CORE().result
        anios = self.CORE().year_list

        img = img.select("LandTrendr")

        # OBTENGO LOS VALORES (ARRAYS)
        ind1 = img.arraySlice(0, 1, 2)
        fit1 = img.arraySlice(0, 2, 3)

        # COMPUTO LA MEDIA DEL INDICE Y LA REPITO PARA QUE SEA UN ARRAY DEL
        # MISMO TAMAÑO
        media = (ind1.arrayReduce(ee.Reducer.mean(), [1])
                 .arrayRepeat(1, ind1.arrayLength(1)))

        # SUMATORIA DE (index - media)**2 Y LO REPITO PARA QUE SEA UN ARRAY
        #  DEL MISMO TAMAÑO
        sstot = (ind1.subtract(media).pow(2)
                 .arrayReduce(ee.Reducer.sum(), [1])
                 .arrayRepeat(1, ind1.arrayLength(1)))

        # SUMATORIA DE (index ajustado - media)**2 Y LO REPITO PARA QUE SEA
        #  UN ARRAY DEL MISMO TAMAÑO
        ssreg = (fit1.subtract(media).pow(2)
                 .arrayReduce(ee.Reducer.sum(), [1])
                 .arrayRepeat(1, ind1.arrayLength(1)))

        # SUMATORIA DE (index - index ajustado)**2 Y LO REPITO PARA QUE
        # SEA UN ARRAY DEL MISMO TAMAÑO
        ssres0 = ind1.subtract(fit1).pow(2).arrayReduce(ee.Reducer.sum(), [1])
        ssres = ssres0.arrayRepeat(1, ind1.arrayLength(1))

        # CREO UN ARRAY CON TANTOS UNOS COMO EL TAMAÑO DEL ARRAY DE INDICES
        arr1 = ee.Image(1).toArray().arrayRepeat(1, ind1.arrayLength(1))

        # DIVIDO
        division = ssres.divide(sstot)

        # A 1 LE RESTO LA DIVISION Y OBTENGO R2
        r2 = arr1.subtract(division)

        # CONVIERTO EN IMAGENES
        imgR2 = r2.arrayReduce(ee.Reducer.first(), [1]).arrayProject(
            [0]).arrayFlatten([["r2"]])

        imgssres = ssres.arrayReduce(ee.Reducer.first(), [1]).arrayProject(
            [0]).arrayFlatten([["ssres"]])

        resultado = namedtuple("Statistics", ["ssres", "r2"])

        return resultado(imgssres, imgR2)

    def slope(self):
        """ Calculate slope of each segment in the LandTrendR fit

        Use:
        LandTrendr(imgcol, maxSegment, index, **kwargs).slope()

        :return: Each image in the collection contains the following bands:

            **0 index**: original index

            **1 {index}_fit**: fitted index (example, ndvi_fit)

            **2 bkp**: breakpoint (1 change, 0 no change)

            **3 year**: image's year

            **4 slope_before**: banda que contiene el slope del tramo anterior al
                año. Si *a1* es el año central, *dif_at = a1-a0*

            **5 slope_after**: banda que contiene el slope del tramo posterior
                al año. Si *a1* es el año central, *dif_ad = a2-a1*

            **6 change**: change's magnitud
                (between -1 (greater loose) and 1 (greater gain))

            **7 angle**: change's angle in radians

        :rtype: ee.ImageCollection

        """

        # APLICO LANDTRENDR
        listaImgsPy = self.breakdown()
        listaImgs = ee.List(listaImgsPy)  # ee.List
        colLTR = ee.ImageCollection(listaImgs)  # ee.ImageCollection

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

    def class1(self, umb_b=0.01, umb_m=0.05):
        """ Método para aplicar una clasificación al resultado de LandTendr

        :Parametros:
        :param umb_b: umbral para determinar el atributo "bajo", el cual
            va entre 0 y este umbral
        :type umb_b: float
        :param umb_m: umbral para determinar los atributos "moderado" y
            "alto", los cuales van entre umb_bajo y umb_mod, y mayor a umb_mod
            respectivamente
        :type umb_m: float

        :return: Una coleccion de imagenes, donde cada imagen tiene una
            banda de nombre "cat", que contiene la categoría, a saber:

            * cat 1: sin grandes cambios (azul)

            * cat 2: perdida suave (naranja)

            * cat 3: perdida alta (rojo)

            * cat 4: ganancia (o recuperacion) suave (amarillo)

            * cat 5: ganancia alta (verde)

            * y tres bandas para la visualizacion
        :rtype: ee.ImageCollection
        """

        col = self.slope()

        def categoria(img):
            d = img.get("system:time_start")
            adelante = ee.String("slope_after")
            atras = ee.String("slope_before")

            umb_bajo = ee.Image.constant(umb_b)
            umb_mod = ee.Image.constant(umb_m)

            at = ee.Image(img).select(atras)
            atAbs = at.abs()

            ad = ee.Image(img).select(adelante)
            adAbs = ad.abs()

            # INTENSIDAD QUIEBRE
            dif = ad.subtract(at)
            media = atAbs.add(adAbs).divide(2)

            # DIRECCION
            at_dec = at.lte(0)  # atras decreciente?
            at_crec = at.gt(0)  # atras creciente?

            ad_dec = ad.lte(0)  # adelante decreciente?
            ad_crec = ad.gt(0)  # adelante creciente?

            # OTRA OPCION DE CATEGORIAS

            # INTENSIDAD

            difAbs = dif.abs()

            menor_bajo = difAbs.lt(umb_bajo)
            menor_mod = difAbs.lt(umb_mod)
            mayor_bajo = difAbs.gte(umb_bajo)
            mayor_mod = difAbs.gte(umb_mod)

            int_baja = menor_bajo
            int_mod = mayor_bajo.And(menor_mod)
            int_alta = mayor_mod

            # int_baja = dif.abs().lt(umb_bajo)
            # int_mod = dif.abs().lt(umb_mod).And(dif.gte(umb_bajo))
            # int_alta = dif.abs().gte(umb_mod)

            cat1 = int_baja  # sin grandes cambios
            cat2 = int_mod.And(dif.lte(0))  # perdida suave
            cat3 = int_alta.And(dif.lte(0))  # perdida alta
            cat4 = int_mod.And(dif.gt(0))  # ganancia suave
            cat5 = int_alta.And(dif.gt(0))  # ganancia alta

            # ESCALO
            cat2 = cat2.multiply(2)
            cat3 = cat3.multiply(3)
            cat4 = cat4.multiply(4)
            cat5 = cat5.multiply(5)

            final = cat1.add(cat2).add(cat3).add(cat4).add(cat5)

            img = img.addBands(
                ee.Image(final).select([0], ["cat"])).addBands(
                ee.Image(at_dec).select([0], ["at_dec"])).addBands(
                ee.Image(at_crec).select([0], ["at_crec"])).addBands(
                ee.Image(ad_dec).select([0], ["ad_dec"])).addBands(
                ee.Image(ad_crec).select([0], ["ad_crec"])).addBands(
                ee.Image(dif).select([0], ["dif"]))

            mask = img.neq(0)

            return img.updateMask(mask).set("system:time_start", d)

        col = col.map(categoria)

        # VISUALIZACION DE IMG
        r = ee.String("viz-red")
        g = ee.String("viz-green")
        b = ee.String("viz-blue")

        def viz(imagen):
            d = imagen.get("system:time_start")
            img = ee.Image(imagen)

            '''
                R     G     B
            1   50,   25,   255   azul
            2   255,  100,  50    naranja
            3   255,  25,   25    rojo
            4   255,  255,  50    amarillo
            5   50,   150,  50    verde
            '''

            cat1 = img.eq(1)
            r1 = cat1.multiply(50)
            g1 = cat1.multiply(25)
            b1 = cat1.multiply(255)

            cat2 = img.eq(2)
            r2 = cat2.multiply(255)
            g2 = cat2.multiply(100)
            b2 = cat2.multiply(50)

            cat3 = img.eq(3)
            r3 = cat3.multiply(255)
            g3 = cat3.multiply(25)
            b3 = cat3.multiply(25)

            cat4 = img.eq(4)
            r4 = cat4.multiply(255)
            g4 = cat4.multiply(255)
            b4 = cat4.multiply(50)

            cat5 = img.eq(5)
            r5 = cat5.multiply(50)
            g5 = cat5.multiply(150)
            b5 = cat5.multiply(50)

            red = r1.add(r2).add(r3).add(r4).add(r5).select([0], [r]).toUint8()

            green = (g1.add(g2).add(g3).add(g4).add(g5)
                     .select([0], [g]).toUint8())

            blue = (b1.add(b2).add(b3).add(b4).add(b5)
                    .select([0], [b]).toUint8())

            return img.addBands(ee.Image([red, green, blue])).set(
                "system:time_start", d)

        col = col.select("cat").map(viz)

        return col

    def classIncendio(self, umbral=0.05, agrupar=False):
        """ Método para aplicar una clasificación al resultado de LandTendr

        :parameter:
        :param umbral: umbra para |pendiente_posterior - pendiente_anterior|
            si pasa este umbral se considera 'cambio'
        :type umbral: float

        :param agrupar: agrupar pixeles
        :type agrupar: bool

        :return: Una coleccion de imagenes, donde cada imagen tiene una banda de
            nombre 'cat', que contiene la categoría, a saber:

            **cat 0**: No incendio

            **cat 1**: posible incendio (perdida y ganancia)

            **cat 2**: posible incendio (solo ganancia)

            tres bandas para la visualizacion y una banda para la intensidad
            denominada "int"
        :rtype: ee.ImageCollection

        """
        col = self.slope()

        def categoria(img):
            d = img.get("system:time_start")

            adelante = ee.String("slope_after")
            atras = ee.String("slope_before")

            indice = ee.Image(img).select(self.index + "_fit")
            umb = ee.Image.constant(umbral)

            # CREO UNA BANDA PARA EL AÑO
            t = ee.Date(ee.Image(img).get("system:time_start")).get("year")
            anio = ee.Image.constant(t).select([0], ["anio"])

            at = ee.Image(img).select(atras)
            # atAbs = at.abs()

            ad = ee.Image(img).select(adelante)
            # adAbs = ad.abs()

            indice_ad = indice.add(ad)
            # indice_at = index.add(at)

            # indice_dif = indice_ad.subtract(index)
            indice_dif = indice_ad.subtract(
                indice)  # .abs().divide(ee.Image.constant(500))

            # INTENSIDAD QUIEBRE
            dif = ad.subtract(at)
            # media = atAbs.add(adAbs).divide(2)

            # DIRECCION
            at_dec = at.lte(0)  # atras decreciente?
            at_crec = at.gt(0)  # atras creciente?

            ad_dec = ad.lte(0)  # adelante decreciente?
            ad_crec = ad.gt(0)  # adelante creciente?

            # OTRA OPCION DE CATEGORIAS

            # INTENSIDAD

            difAbs = dif.abs()

            mayor_mod = difAbs.gte(umb)

            int_alta = mayor_mod

            perdida_alta = int_alta.And(dif.lte(0))  # perdida alta (pa)
            ganancia_alta = int_alta.And(dif.gt(0))  # ganancia alta (ga)

            pa_img = ee.Image(perdida_alta).select([0], ["pa"])
            ga_img = ee.Image(ganancia_alta).select([0], ["ga"])
            ind_dif_img = ee.Image(indice_dif).select([0], ["indice_dif"])

            img = img.addBands(pa_img).addBands(ga_img).addBands(
                ind_dif_img).addBands(anio)

            # mask = img.neq(0)

            return img.set("system:time_start", d)  # .updateMask(mask)

        col = col.map(categoria)

        colL = col.toList(50)
        s = tools.execli(colL.size().getInfo)()
        newCol = ee.List([])

        # ARMO LAS CATEGORIAS TENIENDO EN CUENTA LO QUE SUCEDE ANTES
        # Y DESPUES DEL AÑO EN ANALISIS
        for i in range(0, s):

            ''' ASI ANDA, PRUEBO OTRA COSA Y SINO LO VUELVO A ESTA VER
            img1 = ee.Image(colL.get(i))
            a1 = img1.select("year")
            pa1 = img1.select("pa")
            ga1 = img1.select("ga")
            dif1 = img1.select("indice_dif")	
            
            if (i < s-1):
                img2 = ee.Image(colL.get(i+1))
                a2 = img2.select("year")
                pa2 = img2.select("pa")
                ga2 = img2.select("ga")
                dif2 = img2.select("indice_dif")
            elif (i == s-1):
                # cat 1
                inc = pa1
                cat1 = inc
                
            if (i < s-2):
                img3 = ee.Image(colL.get(i+2))			
                
                a3 = img3.select("year")							
                pa3 = img3.select("pa")
                ga3 = img3.select("ga")
                dif3 = img3.select("indice_dif")
                
                # cat 1
                inc1 = pa1.And(ga2)
                inc2 = pa1.And(ga3)
                inc = inc1.Or(inc2)
                cat1 = inc
                
                # cat 2
                noinc1 = pa1.eq(0) # en la actual no hay perdida
                noinc0 = 
                cat2 = noinc1.And(ga2)
                cat2 = cat2.multiply(ee.Image.constant(2))
                
            elif (i == s-2):
                # cat 1
                inc = pa1.And(ga2)				
                cat1 = inc
                
                # cat 2
                noinc1 = pa1.eq(0)
                cat2 = noinc1.And(ga2)
                cat2 = cat2.multiply(ee.Image.constant(2))
            '''
            img1 = ee.Image(colL.get(i))
            d = img1.get("system:time_start")
            a1 = img1.select("year")
            pa1 = img1.select("pa")
            ga1 = img1.select("ga")
            dif1 = img1.select("indice_dif")

            # EN TODAS LAS IMG EXCEPTO LA PRIMERA, OBTENGO LOS DATOS DE
            # LA IMAGEN ANTERIOR (i-1)
            if i > 0:
                img0 = ee.Image(colL.get(i - 1))
                a0 = img0.select("year")
                pa0 = img0.select("pa")
                ga0 = img0.select("ga")
                dif0 = img0.select("indice_dif")

            # EN TODAS LAS IMG EXCEPTO LA ULTIMA, OBTENGO LOS DATOS DE
            # LA IMAGEN SIGUIENTE (i+1)
            if i < s - 1:
                img2 = ee.Image(colL.get(i + 1))
                a2 = img2.select("year")
                pa2 = img2.select("pa")
                ga2 = img2.select("ga")
                dif2 = img2.select("indice_dif")

            # EN TODAS LAS IMG EXCEPTO LAS ULTIMAS 2, OBTENGO LOS DATOS DE
            # LA IMAGEN SUB SIGUIENTE (i+2)
            if i < s - 2:
                img3 = ee.Image(colL.get(i + 2))

                a3 = img3.select("year")
                pa3 = img3.select("pa")
                ga3 = img3.select("ga")
                dif3 = img3.select("indice_dif")

            # ARMO LAS CATEGORIAS

            # EN LA PRIMER IMAGEN
            if i == 0:
                # 1 GANANCIA EN SIGUIENTE
                cond1 = ga2

                # 2 NADA EN EL SIGUIENTE, GANANCIA EN SUBSIGUIENTE
                # sin_perd_sig = pa2.Not()
                # cond2 = ga3.And(sin_perd_sig)

                cat2 = cond1  # .Or(cond2)

                cat1 = pa1

            # DE LA SEGUNDA A LA ANTEPENULTIMA
            if (i > 0) and (i < s - 2):
                # NADA SIGUIENTE
                sin_gan_sig = ga2.Not()
                sin_perd_sig = pa2.Not()
                nada_sig = sin_gan_sig.And(sin_perd_sig)

                # NADA ANTERIOR
                sin_gan_ant = ga0.Not()
                sin_perd_ant = pa0.Not()
                nada_ant = sin_gan_ant.And(sin_perd_ant)

                # NADA AHORA
                sin_gan = ga1.Not()
                sin_perd = pa1.Not()
                nada = sin_gan.And(sin_perd)

                # CAT 1

                # 1 PERDIDA SEGUIDA DE GANANCIA
                cond1 = pa1.And(ga2)

                # 2 PERDIDA SEGUIDA DE NADA, SEGUIDA DE GANANCIA
                cond2 = pa1.And(nada_sig).And(ga3)

                cat1 = cond1.Or(cond2)

                # CAT 2
                # 1 NADA ANTES, NADA AHORA Y GANANCIA SIGUIENTE
                cat2 = nada_ant.And(nada).And(ga2)

            # EN LA ANTEULTIMA
            if i == s - 2:
                cond1 = pa1.And(ga2)
                cat1 = cond1

                cond3 = pa0.Or(pa1).Not()
                cond4 = ga2.And(cond3)
                cat2 = cond4

            cat2 = cat2.multiply(ee.Image.constant(2))

            final = cat1.add(cat2)
            img1 = img1.addBands(ee.Image(final).select([0], ["cat"]))
            img1 = img1.set("system:time_start", d)
            newCol = newCol.add(img1)

        col = ee.ImageCollection(newCol)

        def viz(img0):
            d = img0.get("system:time_start")

            imagen = ee.Image(img0).select("cat")
            img = ee.Image(imagen)

            indice = ee.Image(img0).select(self.index + "_fit")
            dif = ee.Image(img0).select("indice_dif")
            difA = ee.Image(dif).abs()
            fact = difA.divide(indice)

            # VISUALIZACION DE IMG
            r = ee.String("viz-red")
            g = ee.String("viz-green")
            b = ee.String("viz-blue")
            '''
                R     G     B
            1   255,  25,   25    rojo
            2   255,  100,  50    naranja
            '''

            cat1 = img.eq(1)
            r1 = cat1.multiply(255)
            g1 = cat1.multiply(25)
            b1 = cat1.multiply(25)

            cat2 = img.eq(2)
            r2 = cat2.multiply(255)
            g2 = cat2.multiply(100)
            b2 = cat2.multiply(50)

            red = r1.add(r2).select([0], [r]).toUint8()
            green = g1.add(g2).select([0], [g]).toUint8()
            blue = b1.add(b2).select([0], [b]).toUint8()

            return img.addBands(ee.Image([red, green, blue])).set(
                "system:time_start", d)

        col = col.map(viz)

        if agrupar is False:
            return col
        else:
            def agrupando(img):
                d = ee.Date(img.get("system:time_start"))
                cat = img.select("cat")

                connected = img.toInt().connectedComponents(ee.Kernel.plus(1),
                                                            8).reproject(
                    ee.Projection('EPSG:4326'), None, 30)
                conn2 = connected.mask().select("labels")
                holes = conn2.gt(cat)
                islas = cat.gt(conn2)
                objetos = holes.add(islas)

                kernel = ee.Kernel.plus(1, "pixels", False)
                reducer = ee.Reducer.countDistinct()

                vecinos = (objetos.select(0)
                          .reduceNeighborhood(reducer, kernel, "kernel", False)
                          .reproject(ee.Projection('EPSG:4326'), None, 30))

                objNot = objetos.eq(0)

                final = objNot.Not().reduceNeighborhood(ee.Reducer.max(),
                                                        kernel).reproject(
                    ee.Projection('EPSG:4326'), None, 30)

                final = final.set("system:time_start", d)

                return final

            col = col.map(agrupando)

            return col

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

    def stable_pixels(self, threshold=0.02):
        """ Generate a stable pixels throughout the year image

        :param threshold: value that divides 'change' of 'no change'.
            Defualt: 0.02
        :return: A 8-bits Image with the following bands:

            :cat: pixel's category. Categories are:

            - 1: neutral trend (nor gain neither loss)
            - 2: positive trend (gain)
            - 3: negative trend (loss)

            :viz-red: 8-bits band for category 1
            :viz-green: 8-bits band for category 2
            :viz-blue: 8-bits band for category 3
        """

        col = self.slope()

        # imagen suma de breakpoints
        suma = self.total_bkp(col)

        # sort the collection ascending
        col = col.sort('system:time_start')

        # BANDA indice_fit DE LA PRIMER IMG DE LA COL
        primera = (ee.Image(col.first())
                     .select(self.index + "_fit"))

        ultima = (ee.Image(col.sort('system:time_start', False).first())
                    .select(self.index+'_fit'))

        # OBTENGO LA PENDIENTE GENERAL RESTANDO EL INDICE AL FINAL MENOS
        # EL INDICE INICIAL
        slope_global = (ee.Image(ultima.subtract(primera))
                        .select([0], ["slope"]))

        # CREO LA IMG CON EL UMBRAL
        umb_est = ee.Image.constant(threshold)

        # CREO UNA IMAGEN DE ceros Y CAMBIO EL NOMBRE DE LA BANDA A "bkps"
        ini = ee.Image(0).select([0], ["bkps"])

        # CREO UNA IMG DONDE LOS PIXELES EN LOS QUE
        # |slope_after - slope_before| >= umbral
        # SON unos SINO ceros
        def breakpoints(img, inicial):
            # casteo la imagen inicial
            inicial = ee.Image(inicial)

            # obtengo el valor absoluto del slope
            adelante = img.select("slope_after")  # .abs()
            atras = img.select("slope_before")  # .abs()
            dif = adelante.subtract(atras).abs()

            # creo una imagen binaria
            # si el slope > umbral --> 1
            # sino --> 0
            cond = ee.Image(dif.gte(umb_est))

            return inicial.add(cond)

        bkpImg = ee.Image(col.iterate(breakpoints, ini))

        # estableMask = bkpImg.updateMask(bkpImg.eq(0))

        # DETERMINO LA MASCARA DE PIXELES ESTABLES
        mask = bkpImg.eq(0)

        # ENMASCARO LA IMAGEN DE slopes DEJANDO SOLO LOS ESTABLES
        slope = slope_global.updateMask(mask)

        # slope band to float
        slope = slope.toFloat()

        # DETERMINO UN FACTOR PARA CADA PIXEL
        # DIVIDO POR 0.6 DEBIDO A QUE LA MAXIMA DIFERENCIA
        # ESPERABLE (EN NDVI AL MENOS) ES DE 0.2 (SIN BOSQUE) A
        # 0.8 (CON BOSQUE), LO QUE RESTADO DA 0.6
        slope_fact = slope.abs().divide(0.6)  # .multiply(4)

        # SIN CAMBIOS: PARA LOS QUE LA PENDIENTE DEL SLOPE SEA
        # MENOR AL UMBRAL. QUEDA UNA IMAGEN BINARIA

        estable = (slope.abs()
                   .lte(threshold)
                   .select([0], ["cat"])
                   .toFloat())

        # CRECIMIENTO: PARA LOS QUE LA PENDIENTE ES POSITIVA
        # QUEDA UNA IMG BINARIA MULTIPLICADA X2
        crecimiento = (ee.Image(slope.gte(threshold).multiply(2))
                       .select([0], ["cat"])
                       .toFloat())

        # DECREMENTO: PARA LOS QUE LA PENDIENTE ES NEGATIVA
        # QUEDA UNA IMG BINARIA MULTIPLICADA X3
        decremento = (ee.Image(slope.lt(-threshold).multiply(3))
                      .select([0], ["cat"])
                      .toFloat())

        # IMG CON LAS CATEGORIAS:
        # 'SIN CAMBIO'  = 1
        # 'CRECIMIENTO' = 2
        # 'DECREMENTO'  = 3
        categorias = estable.add(crecimiento).add(decremento)

        # Category band to int 8
        categorias = categorias.toInt8()

        # COLORES
        red = (ee.Image(slope.lt(0))  # IMG CON unos DONDE EL slope < 0
               .multiply(255)  # MULTIPLICA x255
               .multiply(slope_fact)  # MULTIPLICA x slope_fact
               .select([0], ["viz-red"])  # RENOMBRA LA BANDA
               # .toFloat())
               .toInt8())

        green = (ee.Image(slope.gte(0))
                 .multiply(255)
                 .multiply(slope_fact)
                 .select([0], ["viz-green"])
                 # .toFloat())
                 .toInt8())

        blue = (ee.Image(slope.eq(0))
                .multiply(255)
                .select([0], ["viz-blue"])
                # .toFloat())
                .toInt8())

        final = (red.addBands(green)
                 .addBands(blue)
                 .addBands(categorias)
                 .addBands(slope))

        return final

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

    def max_diff(self, category):
        """ Generate an image with the maximum difference found

        :param category: 2 or gain, 3 or loss

        :return: an Image with the following bands:

            :max_change: maximum change value (loss or gain)
            :year: year that occur the max change
        :rtype: ee.Image
        """

        segmentos = self.slope()
        slope_before = segmentos.select("slope_before")
        slope_after = segmentos.select("slope_after")
        change = segmentos.select("change")

        # Determino los pixeles que no contienen un breakpoint para
        # enmascararlos
        # estables, slope = self.stable_pixels()
        # estables = slope.mask().Not()
        estables_mask = self.stable_pixels().mask().Not().select([0])

        # Agrego una banda con el año de la img
        # y enmascaro todas las imgs con la mascara de no stable_pixels
        def agregaFecha(img):
            # mask = img.mask()
            anio = ee.Date(img.get("system:time_start")).get("year")
            imga = ee.Image(anio).select([0], ["year"]).toInt16()
            return img.addBands(imga).updateMask(estables_mask)

        # slope_before = slope_before.map(agregaFecha)
        change = change.map(agregaFecha)  # .select("change")

        # colL = slope_before.toList(50)
        colL = change.toList(50)

        # Primer imagen de la coleccion para la iteracion
        primera = ee.Image(colL.get(0))

        # Resto de la col para la iteracion
        resto = colL.slice(1)

        def maxima(img, first):

            # Casts
            img = ee.Image(img)
            # imagen inicial (primer img de la col)
            firstI = ee.Image(first)

            # test img > img0?

            if category == "loss" or category == 3:
                test = img.select("change").lte(
                    firstI.select("change"))  # binaria
            elif category == "gain" or category == 2:
                test = img.select("change").gte(
                    firstI.select("change"))  # binaria

            # reemplazo condicional

            maxI = firstI.where(test, img)

            # ENMASCARO LA IMAGEN CON LOS PIXELES QUE SE MANTIENEN
            # ESTABLES A LO LARGO DE LA SERIE
            # imgF = maxI.updateMask(stable_pixels)#.addBands(imgy)
            return maxI

        img = ee.Image(resto.iterate(maxima, primera))
        # img = img.updateMask(stable_pixels)
        return img

    def breakpoints(self):
        """ Number of breakpoints

        :return: An object with 2 properties:

            :image: ee.Image containing the following bands:

                :**n_bkp**: number of breakpoints in each pixel

            :total: maximum amount of breakpoints

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

        :return: an Image with the following bands:

            :year_{breakpoint}: year for each breakpoint. Number of breakpoints in
                the whole image may differ from the number of breakpoints in
                each pixel. Therefore, bands will correspond to breakpoint in
                the whole image. Example:
                bkps in the whole image: 5
                bkps in one pixel: 3
                values of bands in that pixel:
                    year_1: 0
                    year_2: 0
                    year_3: 1999
                    year_4: 2005
                    year_5: 2017
            :change_{id}: change value for each breakpoint
            :backward_{id}: backward change value
            :forward_{id}: forward change value
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
        :return: A new class called 'Stretches' with the following properties:

            :img_list: a list of images containing the stretches. Example:
                [img1stretch, img2stretches, img3stretches, img4stretches].
                Each image has the following bands:

                :t{n}_slope: slope of stretch n
                :t{n}_duration: duration of stretch n (in years)
                :t{n}_cat: category for stretch n

            :image: an image with unified results. Will have as many bands as
                found stretches, times 3. In case stretches in one pixel are
                less than stretches in the whole image, the last will be empty.
                Example: The whole image has 4 stretches. The pixel that has 2
                stretches will have the following values:

                :t1_slope: value
                :t2_slope: value
                :t3_slope: 0
                :t4_slope: 0

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

    def classify(self, *args, **kwargs):
        """ Clasificación de los tramos

        :param args: Mismos argumentos que la función stretches()
        :param kwargs: Mismos argumentos que la funcion stretches()
        :return:
        """

        tramos = self.stretches(*args, **kwargs)
        lista = tramos.img_list

        umbral1 = kwargs.get("min_threshold", 0.05)
        umbral2 = kwargs.get("max_threshold", 0.2)

        final = ee.Image.constant(0).select([0],["cat"])

        for n, img in enumerate(lista):
            n += 1

            img = tools.mask2zero(img)

            if n == 1:
                t1_cat = img.select("t1_cat")
                t1_duracion = img.select("t1_duration")
                t1_slope = img.select("t1_slope")

                t1_caida = t1_duracion.multiply(t1_slope).toFloat()
                cat = t1_caida.lte(-umbral1).multiply(3).select([0],["cat"])
                # cat = t1_cat.eq(2).Or(t1_cat.eq(3)).multiply(3).select([0],["cat"])
                final = final.add(cat)
                # continue
                break

            elif n == 2:
                t1_cat = img.select("t1_cat")
                t2_cat = img.select("t2_cat")
                t1_est = t1_cat.eq(1)
                t2_caida = t2_cat.eq(3)
                t2_suave = t2_cat.eq(2)

                cambio = t1_est.And(t2_caida)
                posible = t1_est.And(t2_suave).multiply(2)
                cat = cambio.add(posible).select([0],["cat"])
                final = final.add(cat)
                # continue
                break

            cat0 = ee.Image.constant(0).select([0],["cat"])

            for tramos in range(3, n+1):
                n1 = str(tramos-2)
                n2 = str(tramos-1)
                n3 = str(tramos)

                t1_cat = img.select("t{0}_cat".format(n1))
                t2_cat = img.select("t{0}_cat".format(n2))
                t3_cat = img.select("t{0}_cat".format(n3))

                t2_slope = img.select("t{0}_slope".format(n2))
                t2_duracion = img.select("t{0}_duration".format(n2))

                t1_est = t1_cat.eq(1)
                t1_suave = t1_cat.eq(2)
                t1_caida = t1_cat.eq(3)
                t2_caida = t2_cat.eq(3)
                t2_suave = t2_cat.eq(2)
                t3_est = t3_cat.eq(1)
                t3_crec = t3_cat.eq(4)

                t2_mag = t2_slope.multiply(t2_duracion)

                cambio1 = t1_est.And(t2_caida).And(t3_est)  # stable,
                cambio2 = t1_suave.And(t2_caida).And(t3_est)
                cambio3 = t1_est.And(t2_caida).And(t3_crec)
                cambio4 = t1_suave.And(t2_caida).And(t3_crec)
                cambio5 = t1_caida.And(t2_suave).And(t3_est)
                cambio6 = t1_caida.And(t2_suave).And(t3_crec)
                cambio7 = t1_caida.And(t2_caida).And(t3_est)
                cambio8 = t1_caida.And(t2_caida).And(t3_crec)

                cambio = cambio1.Or(cambio2).Or(cambio3).Or(cambio4).Or(cambio5).Or(cambio6).Or(cambio7).Or(cambio8)

                posible1 = t1_est.And(t2_suave).And(t3_est).multiply(2)
                # posible2 = t1_est.And(t2_mag.lte(-umbral2)).And(t3_est.Or(t3_crec)).multiply(2)

                posible = posible1#.Or(posible2)

                cat = cambio.add(posible).select([0],["cat"])
                cat = tools.mask2zero(cat)

                cat0 = cat0.add(cat)

            final = final.add(cat0)

        return final


'''
if __name__ == "__main__":
    col = ee.ImageCollection("users/rprincipe/AP_tierraDelFuego_mask2number/11nbr_fill")
    lt = LandTrendr(col, "nbr")
    break2 = lt.break2band()
    funciones.asset(break2, "break2band", "users/rprincipe/Pruebas/break2band_1", lt.region)
'''