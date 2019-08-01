import json
from ipywidgets.embed import embed_data
import os, sys
import GenDBScraper.Utilities.nb_utilities as nbu

out_path = '/var/www/sbw25'

if __name__ == "__main__":
    
    strain = sys.argv[1]
    locus_tag = sys.argv[2]

    # Get data from pseudomonas.com and corresponding grid
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
    with open(os.path.join(out_path,'{}_{}.html'.format(strain, locus_tag)), 'w') as fp:
        fp.write(rendered_template)