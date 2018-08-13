# coding=utf-8
""" Module with some tools for IPython Jupyter Notebook and Lab """

from ipywidgets import HTML, Accordion, VBox
from geetools import chart
import ee
ee.Initialize()

def add2map(landtrendr, map):
    """ add a Tab to Map (geetools.ipymap) to plot a LandTrendr result

    :param landtrendr: a landtrendr object
    :type landtrendr: landtrendr.LandTrendr
    """
    region = landtrendr.region
    slope = landtrendr.slope()

    def handler(**kwargs):
        event = kwargs['type']
        coords = kwargs['coordinates']
        wid = kwargs['widget']

        if event == 'click':
            wait = HTML('Loading Chart for point {}..'.format(coords))
            wid.children = [wait]

            region = ee.Geometry.Point(coords)

            ch = chart.Image.series(slope, region,
                                    bands=[landtrendr.index+'_fit',
                                           landtrendr.index],
                                    )

            chart_wid = ch.render_widget()
            wid.children = [chart_wid]

    widget = VBox()

    map.addTab('LandTrendr', handler, widget)


