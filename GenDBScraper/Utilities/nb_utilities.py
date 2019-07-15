
import ipywidgets as widgets
import pandas
import logging
import ipyaggrid
from Bio import SeqIO
from io import StringIO

import re

from GenDBScraper.PseudomonasDotComScraper import PseudomonasDotComScraper, pdc_query

def make_table(strain, locus_tag):
    """ Get the data for strain and locus_tag from pseudomonas.com and render as a table. """

    display(nbu.get_grids(data_tables=nbu.run_pdc(strain, locus_tag)))


def make_table_button(strain, locus_tag):
    """ Return a button. If clicked, display a table for the corresponding data from pdc. """

    def table_button_clicked(b):
        """ Callback for click on the button """
        make_table(strain, locus_tag)

    button = widgets.Button(description=locus_tag)
    button.on_click(table_button_clicked)

    return button

def run_pdc(strain, locus_tag):
    """ Get data for strain and locus tag from pseudomonas.com """


    pdc = PseudomonasDotComScraper(query=pdc_query(strain=strain, feature=locus_tag))
    query_string = "__".join([pdc.query[0].strain, pdc.query[0].feature])
    pdc.connect()
    pdc.run_query()

    results = pdc.results[query_string]

    return results


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

                    pattern = re.compile(r"^Reference$")
                    if pattern.match(column_def['field']):
                        column_def['cellRenderer'] = """function(params) { return '<a http://ncbi.nlm.nih.gov/pubmed/'+params.value+' target=_blank>'+params.value+'</a>'; }"""

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
                        column_def['cellRenderer'] = """function(params) { return '<a href=
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

            js_pre_helpers = []
            js_helpers_custom = """ """
            #js_pre_grid = [""" window.go = gridOptions; window.gd = gridData; function geneCellRenderer(params){ html = '<a href="PA14.ipynb">params.value</a>'; return html; } gridOptions.columnDefs[-1]['cellRenderer']= geneCellRenderer; """]
            js_pre_grid = []
            js_post_grid = []

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
                               js_pre_helpers=js_pre_helpers,
                               js_helpers_custom=js_helpers_custom,
                               js_pre_grid=js_pre_grid,
                               js_post_grid=js_post_grid,
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

# Need to treat each tab and subtabs individually
def get_single_grid(df, title, column_formatting):

    df = df.rename(str, axis='columns')
    grid_options={'columnDefs' : [{'field': c} for c in df.columns],
                          'enableSorting': True,
                          'enableFilter': True,
                          'enableColResize': True,
                          'enableRangeSelection': True,
                          }

    for cd in grid_options['columnDefs']:
        field = cd['field']
        if cd['field'] in column_formatting.keys():
            cd['cellRenderer'] = column_formatting[field]


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

    return grid




def apply_column_formatting(tabs, titles=[], formatting_string=""" """):
    """ Apply the given formatting string to the tab specified by titles

    :param tabs: The tab widget to apply the formatting to.
    :type  tabs: ipywidgets.Tab

    :param titles: Sequence of tab titles and column titles needed to navigate to the tab in question
    :type  titles: list

    :param formatting_string: The formatting string to apply to the specified column.
    :type  formatting_string: str

    """
    t = tabs
    # Navigate to the correct tab by searching for tab titles [vomit].
    for title in titles[:-1]:
        kids = t.children
        # Find index.
        logging.debug("Getting index for title %s", title)
        current_titles = [t.get_title(i) for i in range(len(kids))]
        logging.debug("Current available titles are %s", str(current_titles))
        idx = [ct == title for ct in current_titles].index(True)
        logging.debug("Found idx = %s", str(idx))
        t = kids[idx]

    column_defs = t.grid_options["columnDefs"]

    locate_key = None
    for cd in column_defs:
        if cd['field'] == titles[-1]:
            cd["cellRenderer"] = formatting_string



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


    with open(stdb_results['Network Image'], 'rb') as fp:
        image_widget = widgets.Image(value=fp.read(), format='svg')

    tabs = []
    for key in stdb_results.keys():
        if key == 'Network Image':
            continue
        result = stdb_results[key]
        cds = ColumnDataSource(result)
        data_table = DataTable(source=cds,
                               columns=[TableColumn(field=c, title=c, width=80) for c in list(result.columns)],
                               fit_columns=False
                              )
        tabs.append(Panel(child=data_table, title=key))

    stdb_tabs = Tabs(tabs=tabs)
    display(image_widget)
    show(column(stdb_tabs, width=500))

