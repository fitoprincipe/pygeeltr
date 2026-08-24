#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
from .. import landtrendr
import ee
import pprint
ee.Initialize()

pp = pprint.PrettyPrinter(indent=2)

try:
    from geetools import tools
    exist_geetools = True
except:
    exist_geetools = False


class LandTrendr(unittest.TestCase):

    FOLDER = None
    USER = None

    def output(self, image, name, vis_params):

        folder = LandTrendr.FOLDER if LandTrendr.FOLDER else 'test_landtrendr'

        path = 'users/{}/{}/{}'.format(LandTrendr.USER, folder, name)

        if LandTrendr.USER:
            task = ee.batch.Export.image.toAsset(image, name,
                                                 path,
                                                 region=self.region,
                                                 scale=30)
            task.start()
        else:
            visimage = image.visualize(**vis_params)
            url = visimage.getThumbUrl({'region':self.region})

            print(name, url)

    def output_collection(self, collection, name):
        folder = LandTrendr.FOLDER if LandTrendr.FOLDER else 'test_landtrendr'

        path = 'users/{}/{}/{}'.format(LandTrendr.USER, folder, name)

        if LandTrendr.USER:
            tools.col2asset(collection, path, region=self.region)

    def setUp(self):

        folder = LandTrendr.FOLDER if LandTrendr.FOLDER else 'test_landtrendr'
        # serie = LandTrendr.TIMESERIE if LandTrendr.TIMESERIE else 'test_time_serie'
        serie = 'test_time_serie'

        self.region =  [[[-71.71, -42.77],
                         [-71.71, -42.87],
                         [-71.57, -42.87],
                         [-71.57, -42.77]]]

        self.area = ee.Geometry.Polygon(self.region)  # ee.Geometry.Polygon
        self.center = self.area.centroid()  # ee.Geometry.Point

        time_serie = ee.ImageCollection('users/rprincipe/{}/{}'.format(folder, serie))
        self.principe = landtrendr.LandTrendr.Principe(time_serie, 'nbr', self.area)
        self.kennedy = landtrendr.LandTrendr.Kennedy(time_serie, 'nbr', self.area)
        self.liang = landtrendr.LandTrendr.Liang(time_serie, 'nbr', self.area)

    def test_slope(self):
        slope = self.principe.slopes()
        self.assertEqual(type(slope), ee.ImageCollection)

        self.output_collection(slope, 'test_slope_col')

    def test_statistics(self):
        stats = self.principe.statistics()
        r2 = stats.r2
        ssres = stats.ssres
        self.assertEqual(type(r2), ee.Image)
        self.assertEqual(type(ssres), ee.Image)

        vis_r2 = {'bands':'r2', 'min':0.2, 'max':0.95}

        self.output(r2, 'test_r2', vis_r2)

    def test_breakdown(self):
        bdown = self.principe.breakdown()
        self.assertEqual(type(bdown), list)

        image = bdown[0]
        bands = image.bandNames().getInfo()
        self.assertEqual(bands, ['nbr', 'nbr_fit', 'bkp', 'year'])

        self.output(image, 'test_breakdown', {'bands':['nbr_fit'], 'min':0, 'max':0.8})
        self.output_collection(ee.ImageCollection(ee.List(bdown)), 'test_breakdown_col')

    def test_breakpoints(self):
        breaks = self.principe.breakpoints()
        image = breaks.image
        total = breaks.total

        self.assertEqual(type(image), ee.Image)
        self.assertEqual(type(total), ee.Number)
        self.assertEqual(type(total.getInfo()), int)

        self.output(image, 'test_breakpoints', {'bands':['n_bkp'], 'min':0, 'max':5})

    def test_break2band(self):
        b2b = self.principe.break2band()

        self.assertEqual(type(b2b), ee.Image)
        if exist_geetools:
            values = tools.get_value(b2b, self.center, scale=30, side='client')
            pp.pprint(values)

        self.output(b2b, 'test_break2band', {'bands':['change_0', 'change_1', 'change_2'], 'min':0, 'max':1})

    def test_stretches(self):
        stretches = self.principe.stretches()
        img_list = stretches.img_list
        unified = stretches.image
        one_image = img_list[0]

        self.assertEqual(type(one_image), ee.Image)
        self.assertEqual(type(unified), ee.Image)

        self.output(unified, 'test_stretches', {'bands':['t1_cat'],'min':1,'max':5})
        self.output_collection(ee.ImageCollection(ee.List(img_list)), 'test_stretches_col')

    def test_stack(self):
        stack = self.principe.stack()

        self.assertEqual(type(stack), ee.Image)

        self.output(stack, 'test_stack', {'min':0, 'max':1})

    def test_stable_pixels(self):
        stables = self.principe.stable_pixels()
        self.assertEqual(type(stables), ee.Image)
        self.output(stables, 'test_stables', {'bands':['viz-red',
                    'viz-green', 'viz-blue'], 'min':0, 'max':255})

    def test_total_bkp(self):
        total_bkp = self.principe.total_bkp()
        self.assertEqual(type(total_bkp), ee.Image)
        self.output(total_bkp, 'test_total_bkp', {'bands':['total_bkp'],
                                                  'min':0, 'max':5})

    def test_max_diff_gain(self):
        max_dif = self.principe.max_diff('gain')

        self.assertEqual(type(max_dif), ee.Image)
        self.output(max_dif, 'test_max_gain', {'bands':['max_change'],
                                                  'min':0, 'max':0.5})

    def test_max_diff_loss(self):
        max_dif = self.principe.max_diff('loss')

        self.assertEqual(type(max_dif), ee.Image)
        self.output(max_dif, 'test_max_loss', {'bands':['max_change'],
                                               'min':0, 'max':0.5})


    def test_classify(self):
        classify = self.principe.classify()

        self.assertEqual(type(classify), ee.Image)
        self.output(classify, 'test_classify', {})