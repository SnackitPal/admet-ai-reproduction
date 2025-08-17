
from tdc.single_pred import Tox
import pandas as pd

def download_clintox_data():
    """Downloads the ClinTox dataset from TDC and saves it as a CSV file."""
    data = Tox(name='ClinTox')
    df = data.get_data()
    df.to_csv('admet_predict/data/data.csv', index=False)
    print("ClinTox dataset downloaded and saved to admet_predict/data/data.csv")

if __name__ == "__main__":
    download_clintox_data()
