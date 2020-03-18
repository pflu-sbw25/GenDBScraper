import ipywidgets as widgets
import pandas
import logging
import ipyaggrid
from Bio import SeqIO
from io import StringIO

import re
import json

from IPython.display import IFrame, clear_output, Image
from GenDBScraper.PseudomonasDotComScraper import PseudomonasDotComScraper, pdc_query

from GenDBScraper.StringDBScraper import StringDBScraper, stringdb_query


def sbw25_okm():
    clear_output()
    frame = IFrame("https://openknowledgemaps.org/map/4aafb7d70516de0f56190d374bf398c8&embed=true", width=1000, height=1000)
    return frame
    #display(frame)
    
def feature_okm_js(locus_tag):
    okms = pandas.read_json("pflu_okm_urls_20190424.json", typ='series', orient='records')

    if locus_tag.upper() in okms.keys():
        frame = '''<iframe width="1200" height="720" src="{0:s}"></iframe>'''.format(okms[locus_tag.upper()]+'''&embed=true''')

    else:
        frame =  "<iframe></iframe>"

    return frame

def feature_okm(locus_tag):
    okms = pandas.read_json("pflu_okm_urls_20190424.json", typ='series', orient='records')

    if locus_tag.upper() in okms.keys():
        frame = IFrame(okms[locus_tag.upper()]+"&embed=true", width=900, height=800)

    else:
        frame = widgets.Label("No maps found for {}".format(locus_tag))

    clear_output(wait=True)    
    #display(frame)
    return frame

def run_pdc(strain=None, locus_tag=None):
    
    clear_output(wait=True)
    query=pdc_query(strain=strain, feature=locus_tag)
    pdc = PseudomonasDotComScraper(query)
    pdc.connect()
    pdc.run_query()
    
    
    results = pdc.results[query_string]
    
    return results
    

# Need to treat each tab and subtabs individually
def get_grids(data_tables):
    """ Create grid view of all data tables"""
    
    if not isinstance(data_tables, dict):
        raise TypeError("Input parameter 'data_tables' must be of type dict. Received type is {}".format(type(data_tables)))
        
    tabs = widgets.Tab() 
    
    children = []
    titles = []
    
    skipped = ["Ortholog xml"]
    
    for i, title in enumerate(data_tables.keys()):
        if title in skipped:
            logging.debug("Skipping %s", title)
            continue
            
        df = data_tables[title]
        if df is None:
            logging.debug("Skipping %s", title)
            continue

        if isinstance(df, pandas.DataFrame):    
            if df.empty:
                logging.debug("Skipping %s", title)
                continue
                
            df = df.rename(str, axis='columns')
            
            grid_options={'columnDefs' : [{'field': c} for c in df.columns],
                          'enableSorting': True,
                          'enableFilter': True,
                          'enableColResize': True,
                          'enableRangeSelection': True,
                          }
          
            if title.lower() == "ortholog group":                           
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^GI$")            
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://www.ncbi.nlm.nih.gov/protein/'+params.value+' target=_blank>'+params.value+'</a>'; }"""
          
            if title.lower() == "ortholog cluster":                           
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^GI \(Strain [1,2]\)$")            
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://www.ncbi.nlm.nih.gov/protein/'+params.value+' target=_blank>'+params.value+'</a>'; }"""
          
            if title.lower() == "cross-references":
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^[U,u]rl$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href='+params.value+' target=_blank>'+params.value+'</a>'; }"""
            
            if title.lower() == "individual mappings":
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^PMID$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://ncbi.nlm.nih.gov/pubmed/'+params.value+' target=_blank>'+params.value+'</a>'; }"""
            
            if title.lower() == "gene ontology":
                for column_def in grid_options['columnDefs']:
                    # GO Accession
                    pattern = re.compile(r"^Accession$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://www.ebi.ac.uk/QuickGO/GTerm?id='+params.value+' target=_blank>'+params.value+'</a>'; }"""
                    
                    # ECO accession
                    pattern = re.compile(r"^Evidence Ontology \(ECO\) Code$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://www.ebi.ac.uk/ontology-lookup/?termId='+params.value+' target=_blank>'+params.value+'</a>'; }"""
    
            if title.lower() == "functional predictions from interpro":
                for column_def in grid_options['columnDefs']:
                   
                    pattern = re.compile(r"^Interpro Accession$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://www.ebi.ac.uk/interpro/entry/'+params.value+' target=_blank>'+params.value+'</a>'; }"""
           
            if re.match(r'^transposon.*$', title.lower() ):
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^Reference$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://ncbi.nlm.nih.gov/pubmed/'+params.value+' target=_blank>'+params.value+'</a>'; }"""

            if title.lower() == 'genes':
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^Unnamed: 7$",flags=re.IGNORECASE)
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { 
    let v = params.value;
    function clicked(){
        let new_cell = Jupyter.notebook.insert_cell_below().set_text("This feature is not implemented yet.");
    }

    let b = document.createElement('button');
    b.innerHTML = v;
    b.style = "background-color:bisque; margin:1px 10px 1px 2px;";
    b.title = "Open gene table";
    b.addEventListener("click", function (){clicked()}, false);
    // b.addEventListener("click", function (){clicked()}, false);

    return b;
} """
            
            if title.lower() == 'references':
                for column_def in grid_options['columnDefs']:
                    pattern = re.compile(r"^Pubmed_id$",flags=re.IGNORECASE)
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a href=http://ncbi.nlm.nih.gov/pubmed/'+params.value+' target=_blank>'+params.value+'</a>'; }"""
            
            
            g = ipyaggrid.Grid(grid_data = df,
                               grid_options=grid_options,
                               center=False,
                               theme='ag-theme-fresh',
                               grid_options_multi=[],
                               columns_fit='',
                               index=True,
                               keep_multiindex=False,
                               compress_data=True,
                               quick_filter=True,
                               export_csv=True,
                               export_excel=True,
                               show_toggle_delete=False,
                               show_toggle_edit=False,
                               paste_from_excel=True,
                               export_mode='disabled',
                               export_to_df=True,
                               hide_grid=False,
                               menu=None,
                               )            
            
            children.append(g)
            
        elif isinstance(df, dict):
            if df == {}:
                logging.debug("Skipping %s", title)
            
            g = get_grids(df)
            children.append(g)              
            
        elif isinstance(df, list):
            if len(df) == 0:
                logging.debug("Skipping %s", title)
            continue         
        
        titles.append(title) 
            
    tabs.children = children
    
    assert len(children) == len(titles)
    
    for i, title in enumerate(titles):
        tabs.set_title(i, title)

    return tabs

def run_stdb(locus_tag):
    clear_output(wait=True)

    gene_sub_pattern = re.compile(r'([a-z](?=[0-9]))')
    gene=gene_sub_pattern.sub(r'\1_', locus_tag)

    stdb = StringDBScraper(query=stringdb_query(taxonId=216595, features=[gene]))

    stdb.connect()

    stdb.update_features()

    stdb_results = dict()

    stdb_results['Network Image'] = stdb.network_image()
    stdb_results['Network Interactions'] = stdb.network_interactions()
    stdb_results['Interaction Partners'] = stdb.interaction_partners(required_score=300)
    stdb_results['Functional Enrichments'] = stdb.functional_enrichments()
    stdb_results['Interaction Enrichments'] = stdb.interaction_enrichments()

    return stdb_results

def get_stdb_grids(stdb_results):
    
    tabs = widgets.Tab() 
    children = []
    titles = []
    
    skipped = []
    
    for i, key in enumerate(stdb_results.keys()):
        if key in skipped:
            logging.debug("Skipping %s", key)
            continue
            
        if key == 'Network Image':
            with open(stdb_results['Network Image'], 'rb') as fp:
                image_widget = widgets.Image(value=fp.read(), format='svg')

            children.append(image_widget)
            titles.append(key)
            continue
            
        df = stdb_results[key]
        
        grid_options={'columnDefs' : [{'field': c} for c in df.columns],
                          'enableSorting': True,
                          'enableFilter': True,
                          'enableColResize': True,
                          'enableRangeSelection': True,
                          }

        grid =  ipyaggrid.Grid(grid_data = df,
                               grid_options=grid_options,
                               center=False,
                               theme='ag-theme-fresh',
                               grid_options_multi=[],
                               columns_fit='',
                               index=False,
                               keep_multiindex=False,
                               compress_data=True,
                               quick_filter=True,
                               export_csv=True,
                               export_excel=True,
                               show_toggle_delete=False,
                               show_toggle_edit=False,
                               paste_from_excel=True,
                               export_mode='disabled',
                               export_to_df=True,
                               hide_grid=False,
                               menu=None,
                               )
        children.append(grid)
            
        
        titles.append(key) 
            
    tabs.children = children
    
    assert len(children) == len(titles)
    
    for i, title in enumerate(titles):
        tabs.set_title(i, title)
                  
    return tabs

    stdb_tabs = Tabs(tabs=tabs)
    
    return stdb_tabs

def html_template(strain, locus_tag):
    
    ht = """
<html>
  <head>

    <title>Feature results</title>

    <!-- Load RequireJS, used by the IPywidgets for dependency management -->
    <script 
      src="https://cdnjs.cloudflare.com/ajax/libs/require.js/2.3.4/require.min.js" 
      integrity="sha256-Ae2Vz/4ePdIu6ZyI/5ZGsYnb+m0JlOmKPjt6XZ9JJkA=" 
      crossorigin="anonymous">
    </script>

    <!-- Load IPywidgets bundle for embedding. -->
    <script
      data-jupyter-widgets-cdn="https://cdn.jsdelivr.net/npm/"
      src="https://unpkg.com/@jupyter-widgets/html-manager@*/dist/embed-amd.js" 
      crossorigin="anonymous">
    </script>

    <!-- The state of all the widget models on the page -->
    <script type="application/vnd.jupyter.widget-state+json">
      {manager_state}
    </script>
  </head>

  <body>

    <h2>Data from pseudomonas.com</h2>

    <div id="pdc-widget">
      <!-- This script tag will be replaced by the view's DOM tree -->
      <script type="application/vnd.jupyter.widget-view+json">
        {widget_views[0]}
      </script>
    </div>

    <hrule />
    <h2>Data from strings-db</h2>
    <div id="stdb-widget">
      <!-- This script tag will be replaced by the view's DOM tree -->
      <script type="application/vnd.jupyter.widget-view+json">
        {widget_views[1]}
      </script>
    </div>
    
    <hrule />
    <h2> Open Knowledge Map </h2>
    <div id="okm">
    {okm}
    </div>
  </body>
</html>
"""
    
    return ht
