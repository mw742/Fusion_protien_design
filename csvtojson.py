import csv
import json
import pandas as pd

csv_data=pd.read_csv("/users/wmm/desktop/result.csv", sep = ",", low_memory=False)
csv_data.to_json("test.json", orient = "records")