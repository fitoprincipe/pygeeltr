""" Helpers """
import ee


def filter_magnitude(collection, magnitude, index='nbr'):
    """ Mask out pixels with magnitude greater or equal than first item of
    tuple and lower than the last item """
    band = '{}_magnitude'.format(index)
    values = ee.List(list(magnitude))
    def wrap(image):
        mag = image.select(band)
        mask = mag.gte(ee.Number(values.get(0))).And(mag.lt(ee.Number(values.get(1))))
        return image.updateMask(mask)
    return collection.map(wrap)


def filter_duration(collection, duration, index='nbr'):
    """ Mask out pixels with duration greater or equal than first item of
    tuple and lower than the last item """
    band = '{}_duration'.format(index)
    values = ee.List(list(duration))
    def wrap(image):
        mag = image.select(band)
        mask = mag.gte(ee.Number(values.get(0))).And(mag.lt(ee.Number(values.get(1))))
        return image.updateMask(mask)
    return collection.map(wrap)