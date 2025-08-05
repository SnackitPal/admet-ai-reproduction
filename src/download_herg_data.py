
from tdc.single_pred import Tox
import pandas as pd

def download_herg_data():
    """Downloads the hERG dataset from TDC and saves it as a CSV file."""
    data = Tox(name='hERG')
    df = data.get_data()
    df.to_csv('data/herg_data.csv', index=False)
    print("hERG dataset downloaded and saved to data/herg_data.csv")

if __name__ == "__main__":
    download_herg_data()
