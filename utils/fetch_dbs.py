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

import gzip
from collections import defaultdict
from io import BytesIO

import pandas as pd
import requests


def get_geo_family_soft_samples(
    url: str,
    trim_prefix: bool = True,
    collapse_repeated: bool = True,
    repeated_sep: str = " | ",
) -> pd.DataFrame:
    """
    Download a GEO family.soft.gz file from a URL and parse sample metadata
    into a DataFrame with one row per sample.

    Parameters
    ----------
    url : str
        URL to the GEO family.soft.gz file.
    trim_prefix : bool, default True
        Remove '!Sample_' prefix from column names.
    collapse_repeated : bool, default True
        If a metadata field appears multiple times for the same sample
        (e.g. characteristics_ch1), join values into one string.
        If False, repeated fields get numbered suffixes.
    repeated_sep : str, default " | "
        Separator used when collapsing repeated fields.

    Returns
    -------
    pd.DataFrame
        One row per sample, columns are sample metadata fields.
    """
    response = requests.get(url)
    response.raise_for_status()

    with gzip.GzipFile(fileobj=BytesIO(response.content)) as gz:
        samples = []
        current = None
        in_sample_table = False

        for raw_line in gz:
            line = raw_line.decode("utf-8", errors="replace").rstrip("\n")

            # Start of a sample block
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                sample_id = line.split(" = ", 1)[1].strip()
                current = {"sample_id": sample_id}
                in_sample_table = False
                continue

            # Ignore everything until the first sample block
            if current is None:
                continue

            # Skip sample expression table contents
            if line.startswith("!sample_table_begin"):
                in_sample_table = True
                continue
            if line.startswith("!sample_table_end"):
                in_sample_table = False
                continue
            if in_sample_table:
                continue

            # Sample metadata lines
            if line.startswith("!Sample_"):
                key, value = (line.split(" = ", 1) + [""])[:2]

                if trim_prefix:
                    key = key.removeprefix("!Sample_")

                value = value.strip()

                if collapse_repeated:
                    if key in current and current[key]:
                        current[key] = f"{current[key]}{repeated_sep}{value}"
                    else:
                        current[key] = value
                else:
                    if key not in current:
                        current[key] = value
                    else:
                        i = 2
                        while f"{key}_{i}" in current:
                            i += 1
                        current[f"{key}_{i}"] = value

        if current is not None:
            samples.append(current)

    return pd.DataFrame(samples)