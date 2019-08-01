
from ipywidgets.embed import embed_data
import os, sys
import logging
import json
import GenDBScraper.Utilities.nb_utilities as nbu
from multiprocessing import Pool

OUT_PATH = '/var/www/sbw25'

def process_tag(tag):
    strain = "sbw25"
    locus_tag = r'pflu{0:04d}'.format(tag)

    logging.info("Processing %s .", locus_tag)

    # Get data from pseudomonas.com and corresponding grid
    try:
        pdc_grid = nbu.get_grids(nbu.run_pdc(strain, locus_tag))

        # Get data from strings-db and the corresponding grid
        stdb_grid = nbu.get_stdb_grids(nbu.run_stdb(locus_tag))

        #okm = nbu.feature_okm(locus_tag)
        #data = embed_data(views=[pdc_grid, stdb_grid, okm])

        data = embed_data(views=[pdc_grid, stdb_grid])


        manager_state = json.dumps(data['manager_state'])
        widget_views = [json.dumps(view) for view in data['view_specs']]

        rendered_template = nbu.html_template(strain, locus_tag).format(manager_state=manager_state,
                                                                        widget_views=widget_views,
                                                                        okm = nbu.feature_okm_js(locus_tag),
                                                                        )
        with open(os.path.join(OUT_PATH,'{}_{}.html'.format(strain, locus_tag)), 'w') as fp:
            fp.write(rendered_template)
        return 0

    except:
        return 1

if __name__ == "__main__":

    #tags = range(1,6102)
    tags = range(3957, 6102)

    nproc = 20

    pool = Pool(nproc)

    print(pool.map(process_tag, tags))

