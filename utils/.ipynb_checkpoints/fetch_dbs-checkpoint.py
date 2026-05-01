import pandas as pd
import gzip
import requests
from io import BytesIO

def get_EGAD_files(EGAD: str)-> pd.DataFrame:
    response = requests.get(f"https://metadata.ega-archive.org/datasets/{EGAD}/files")
    EGAFs = response.json()

    data = []
    for EGAF in EGAFs:
        #print(EGAF)
        EGANs = requests.get(f"https://metadata.ega-archive.org/files/{EGAF['accession_id']}/samples").json()
        EGAXs = requests.get(f"https://metadata.ega-archive.org/files/{EGAF['accession_id']}/experiments").json()
        
        # only look on files assotiated with experiments
        if len(EGAXs) == 0:
            continue
            
        if len(EGANs) > 1:
            raise Exception(f"File ({EGAF['accession_id']} has multiple EGANs ({[n['accession_id'] for n in EGANs]}) only first one will be used")
        if len(EGAXs) > 1:
            raise Exception(f"File ({EGAF['accession_id']} has multiple EGAXs ({[x['accession_id'] for x in EGAXs]}) only first one will be used")
        
        
        data.append([
            EGAD,
            EGAF['accession_id'],                     #EGAF
            EGAXs[0]['accession_id'],                 #EGAX
            EGANs[0]['accession_id'],                 #EGAN
            EGANs[0]["title"],                        #AUTHOR_ID
            EGAXs[0]['library_construction_protocol'],#LIBRARY_TYPE
            EGANs[0]['description'],                   
        ])
    df = pd.DataFrame(data, columns=["EGAD","EGAF","EGAX","EGAN","AUTHOR_ID","LIBRARY_TYPE",'DESCRIPTION'])
    return df

def get_geo_sample_metadata(url: str) -> pd.DataFrame:
    """
    Download a gzipped GEO series matrix file from a URL,
    extract lines starting with '!Sample_', and return as a DataFrame.

    Parameters
    ----------
    url : str
        URL to the .gz GEO series matrix file.
    as_header : bool, default True
        If True, use the first row as column headers.
    trim_prefix : bool, default True
        If True, remove '!Sample_' prefix from column names.

    Returns
    -------
    pd.DataFrame
    """
    response = requests.get(url)
    response.raise_for_status()

    with gzip.GzipFile(fileobj=BytesIO(response.content)) as gz:
        lines = [
            line.decode("utf-8").strip()
            for line in gz
            if line.decode("utf-8").startswith("!Sample_")
        ]

    data = [line.split("\t") for line in lines]
    df = pd.DataFrame(data)
    df = df.T
    
    df.columns = df.iloc[0]
    df = df[1:].reset_index(drop=True)
    df.columns = [
              col.replace("!Sample_", "") if isinstance(col, str) else col
              for col in df.columns
    ]

    return df