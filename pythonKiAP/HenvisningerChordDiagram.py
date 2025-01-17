import pandas as pd
import numpy as np
from d3blocks import D3Blocks  
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
import os
import HenvisningerChordDiagram as hcd

#--------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    print(" - Loading referreal data from CSV")
    filetype = 'fullcsv'
    if filetype == 'Short':
        input_file = "C:/Users/kbschmidt/OneDrive - KiAP Fonden/Documents/SQL/Henvisninger_SORberiget_short.txt"
    elif filetype == 'fullcsv':
        input_file = "C:/Users/kbschmidt/OneDrive - KiAP Fonden/Documents/SQL/Henvisninger_SORberiget.txt"
        #rember to add column headers manually when saving from SSMS

    df_ref = pd.read_csv(input_file, sep='\t') 

    Nrows  = len(df_ref['aar'])

    #cols_to_show = 'HealthInstitutionEntityTypeName'
    #cols_to_show = 'InstitutionOwnerEntityTypeName'
    cols_to_show = 'EntityName'

    piv_source_col = 'afs_'+cols_to_show
    piv_target_col = 'mod_'+cols_to_show
    piv_values_col = 'antal'

    print(" - Pivoting input data ")
    df_piv = pd.pivot_table(df_ref, index=piv_source_col, columns=piv_target_col, values=piv_values_col, aggfunc='sum').fillna(0)

    print(" - Converting pivot table to format required by Chord diagram function")
    # Convert pivot table to long format
    df_long = df_piv.reset_index().melt(id_vars=piv_source_col, var_name=piv_target_col, value_name=piv_values_col)
    # Rename columns to match D3Blocks requirements
    df_long.rename(columns={piv_source_col: 'source', piv_target_col: 'target', piv_values_col: 'weight'}, inplace=True)
    # Normalize values
    #df_long['weight'] = df_long['weight'] / df_long['weight'].max()
    # Filter out zero or irrelevant values (optional)
    df_long = df_long[df_long['weight'] > 0]

    print(" - Initialize and plot Chord diagram ")
    d3 = D3Blocks(chart='Chord', frame=False)

    d3.set_node_properties(df_long, opacity=0.7, cmap='tab20')
    d3.set_edge_properties(df_long, color='source', opacity='source')

    #d3.node_properties.get('kommune')['color']='#ff0000'
    #d3.node_properties.get('kommune')['opacity']=1

    # Make edits to highlight the Nuclear Edge
    d3.edge_properties.loc[(d3.edge_properties['source'] == 'kommune') & (d3.edge_properties['target'] == 'privat'), 'color'] = '#ff0000'

    # Show the chart
    htmlpathname = 'C:/Users/kbschmidt/OneDrive - KiAP Fonden/Pictures/PythonFigures/henvisning_chord_'+cols_to_show+'.html'
    d3.show(showfig=False, filepath=htmlpathname, title='Chord of col:'+cols_to_show,save_button=False)

    htmltext ='Kolonner for modtager og afsender af henvisning: {0}<br> Antal rækker i input henvisningsdata: {1}<br> Antal rækker i dataframe når henvisningspar med antal=0 fjernes: {1}<br>'.format(cols_to_show,str(Nrows),str(len(df_long)))
    hcd.add_html_textbox(htmlpathname,htmltext,verbose=True)

    #hcd.save_using_selenium(htmlpathname)

#--------------------------------------------------------------------------------------------------------------------
def add_html_textbox(htmlpathname,htmltext,verbose=True):
    """
    Add a textbox to an html file by replace "body"
    And remove add-box in chord
    """
    if verbose: print(' - Loading HTML file to add text box with info ')
    # Open the html file to update
    with open(htmlpathname, "r", encoding="UTF-8") as file:
        html_content = file.read()

    # Add the info box
    info_box = """
    <div id="info-box" style="position: absolute; bottom: 20px; left: 20px; background: #f9f9f9; padding: 10px; border: 1px solid #ccc; border-radius: 5px;">
        <p><strong>Diagram Info:</strong></p>
        <p>{0}</p>
    </div>
    """.format(htmltext)
    html_content = html_content.replace("</body>", f"{info_box}\n</body>")

    if verbose: print(' - Removing add-box')
    html_content = html_content.replace("<script async src='https://media.ethicalads.io/media/client/ethicalads.min.js'></script>","")
    html_content = html_content.replace("<div data-ea-publisher='erdogantgithubio' data-ea-type='text' data-ea-style='stickybox'></div>","")

    if verbose: print(' - Save the updated HTML file')
    with open(htmlpathname, "w", encoding="UTF-8") as file:
        file.write(html_content)

#--------------------------------------------------------------------------------------------------------------------
def save_using_selenium(htmlpathname):
    """
    Using Selenium to excecute the javascript code piece, that saves the figure from the HTML file
    """
    # Path to your ChromeDriver (update this to your path)
    chromedriver_path = "C:/Program Files/Google/Chrome/Application/chrome.exe"

    # Configure Chrome to run in headless mode
    options = Options()
    options.add_argument("--headless")  # Run in headless mode
    options.add_argument("--disable-gpu")  # Disable GPU rendering

    # Set up ChromeDriver
    service = Service(chromedriver_path)
    driver = webdriver.Chrome(service=service, options=options)

    # Path to the HTML file
    html_file_path = os.path.abspath(htmlpathname)

    # Load the HTML file in the browser
    driver.get(f"file://{html_file_path}")

    # Execute the JavaScript code to save the SVG
    js_code = """
        var svgData = document.querySelector('svg').outerHTML;
        var blob = new Blob([svgData], {type: 'image/svg+xml;charset=utf-8'});
        var fileReader = new FileReader();
        fileReader.onload = function(event) {
            var svgContent = event.target.result;
            return svgContent;
        };
        fileReader.readAsText(blob);
    """
    svg_content = driver.execute_script(js_code)

    # Save the SVG content to a file
    with open(htmlpathname.replace('html','svg'), "w", encoding="utf-8") as svg_file:
        svg_file.write(svg_content)

    # Close the browser
    driver.quit()
#--------------------------------------------------------------------------------------------------------------------