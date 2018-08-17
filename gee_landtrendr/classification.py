# coding=utf-8
""" Module holding methods to classify landtrendr results """

import ee
ee.Initialize()

from geetools import tools
from collections import namedtuple


def classify(ltr, *args, **kwargs):
    """ Clasificación de los tramos

    :param args: Mismos argumentos que la función stretches()
    :param kwargs: Mismos argumentos que la funcion stretches()
    :return:
    """

    tramos = stretches(*args, **kwargs)
    lista = tramos.img_list

    umbral1 = kwargs.get("min_threshold", 0.05)
    umbral2 = kwargs.get("max_threshold", 0.2)

    final = ee.Image.constant(0).select([0],["cat"])

    for n, img in enumerate(lista):
        n += 1

        # img = tools.mask2zero(img)
        img = img.unmask()

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


def class1(ltr, umb_b=0.01, umb_m=0.05):
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

    col = ltr.slope

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


def classIncendio(ltr, umbral=0.05, agrupar=False):
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
    col = ltr.slope()

    def categoria(img):
        d = img.get("system:time_start")

        adelante = ee.String("slope_after")
        atras = ee.String("slope_before")

        indice = ee.Image(img).select(ltr.fit_band + "_fit")
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

        indice = ee.Image(img0).select(ltr.fit_band + "_fit")
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


def classIncendioImage(classified_image, time_list):
    """ Create a unique image with encoded values for the fire year of
    occurrence """
    pass


def stable_pixels(ltr, threshold=0.02):
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

    col = ltr.slope()

    # imagen suma de breakpoints
    suma = ltr.total_bkp(col)

    # sort the collection ascending
    col = col.sort('system:time_start')

    # BANDA indice_fit DE LA PRIMER IMG DE LA COL
    primera = (ee.Image(col.first())
               .select(ltr.fit_band + "_fit"))

    ultima = (ee.Image(col.sort('system:time_start', False).first())
              .select(ltr.fit_band + '_fit'))

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


def max_diff(ltr, category):
    """ Generate an image with the maximum difference found

    :param category: 2 or gain, 3 or loss

    :return: an Image with the following bands:

        :max_change: maximum change value (loss or gain)
        :year: year that occur the max change
    :rtype: ee.Image
    """

    segmentos = ltr.slope()
    slope_before = segmentos.select("slope_before")
    slope_after = segmentos.select("slope_after")
    change = segmentos.select("change")

    # Determino los pixeles que no contienen un breakpoint para
    # enmascararlos
    # estables, slope = ltr.stable_pixels()
    # estables = slope.mask().Not()
    estables_mask = ltr.stable_pixels().mask().Not().select([0])

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

def stretches(ltr, min_threshold=0.05, max_threshold=0.2):
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
    bkps = ltr.breakpoints
    img_bkp = bkps.image

    # Total number of breakpoints
    total_bkp = bkps.total.getInfo()
    max_tramos = total_bkp - 1

    bandas_a = ["year_"+str(i) for i in range(total_bkp)]
    bandas_fit = ["fit_"+str(i) for i in range(total_bkp)]
    bandas = bandas_a + bandas_fit

    b2b = ltr.break2band().select(bandas)

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
            # funciones.asset(img, nombre, "users/rprincipe/Pruebas/"+nombre, ltr.region)
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

        # funciones.asset(img, nombre, "users/rprincipe/Pruebas/"+nombre, ltr.region)

        # bandas = img.getInfo()["bands"]
        # print "tramos:", tramos, [i["id"] for i in bandas]

        img_final = img_final.add(img)

    # return

    resultado = namedtuple("Stretches", ["img_list", "image"])

    return resultado(listaimg, img_final)