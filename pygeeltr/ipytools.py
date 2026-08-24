# coding=utf-8
""" Module with some tools for IPython Jupyter Notebook and Lab """

import ee


def add2map(landtrendr, map, scale=30, name='{band}_ltr', range=None):
    """ add a Tab to Map (geetools.ipymap) to plot a LandTrendr result

    :param landtrendr: a landtrendr object
    :type landtrendr: landtrendr.LandTrendr
    """
    from ipywidgets import HTML, VBox
    import ipygee as ui

    # TODO: make it async so it does not block other widgets
    bands = [landtrendr.fit_band+'_fit', landtrendr.fit_band]
    col = landtrendr.breakdown  # .select(bands)

    def handler(**kwargs):
        event = kwargs['type']
        coords = kwargs['coordinates']
        wid = kwargs['widget']

        if event == 'click':
            wait = HTML('Loading Chart for point {}..'.format(coords))
            wid.children = [wait]

            region = ee.Geometry.Point(coords)

            try:
                params = dict(imageCollection=col,
                              region=region,
                              bands=bands,
                              scale=scale)
                if range:
                    params['range'] = range

                ch = ui.chart.Image.series(**params)
                ch.title = 'LandTrendr fit\n for point {}'.format(coords)
                chart_wid = ch.renderWidget()
            except Exception as e:
                chart_wid = HTML('{}'.format(e))

            wid.children = [chart_wid]

    widget = VBox()
    name = name.format(band=landtrendr.fit_band)
    map.addTab(name, handler, widget)


