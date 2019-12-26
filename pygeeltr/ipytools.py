# coding=utf-8
""" Module with some tools for IPython Jupyter Notebook and Lab """

from ipywidgets import HTML, VBox
from geetools.ui import chart
import ee
ee.Initialize()


def add2map(landtrendr, map):
    """ add a Tab to Map (geetools.ipymap) to plot a LandTrendr result

    :param landtrendr: a landtrendr object
    :type landtrendr: landtrendr.LandTrendr
    """
    # TODO: make it async so it does not block other widgets
    bands = [landtrendr.fit_band+'_fit', landtrendr.fit_band]
    slope = landtrendr.slopes  # .select(bands)

    def handler(**kwargs):
        event = kwargs['type']
        coords = kwargs['coordinates']
        wid = kwargs['widget']

        if event == 'click':
            wait = HTML('Loading Chart for point {}..'.format(coords))
            wid.children = [wait]

            region = ee.Geometry.Point(coords)

            try:
                ch = chart.Image.series(slope, region,
                                        bands=bands)
                ch.title = 'LandTrendr fit\n for point {}'.format(coords)
                chart_wid = ch.render_widget()
            except Exception as e:
                chart_wid = HTML('{}'.format(e))

            wid.children = [chart_wid]

    widget = VBox()

    map.addTab('LandTrendr', handler, widget)


