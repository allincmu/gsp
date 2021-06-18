import pandas as pd

def add_xls(to_xls, panel_name, num_panels, deliv_date):
  to_xls["panel name"].append(panel_name)
  to_xls["num panels"].append(num_panels)
  to_xls["deliv date"].append(deliv_date)

df = pd.read_excel('Copy of EnclosurePriorityList 05.17.21 (002).xlsx', index_col = 0)

buildings = df.to_dict(orient='dict')
sbc = buildings["NEW SBC Enclosure"]
TC = buildings["Telecom"]
ICC = buildings["ICC"]
HHW = buildings["HHW"]
SD = buildings["Start Date"]

SBC_Tot = 0
TC_Tot = 0
ICC_Tot = 0
HHW_Tot = 0

to_xls = {"panel name": [], "num panels" : [], "deliv date" : []}

for building in sbc:
  if sbc[building] != 0 and sbc[building] != "NONE" and sbc[building] != "UPSIZE on site":
    panel_name = (str(building) + " - New SBC Enclosure")
    add_xls(to_xls, panel_name, sbc[building], SD[building])
    SBC_Tot += 1
  else:
    print("\trejected: " + str(sbc[building]))


  if TC[building] != 0 and  TC[building] != "NONE":
    panel_name = (str(building) + " - Telecom Enclosure")
    add_xls(to_xls, panel_name, TC[building], SD[building])
    TC_Tot += 1
  else:
    print("\trejected: " + str(TC[building]))


  if ICC[building] != 0 and ICC[building] != "NONE":
    panel_name = (str(building) + " - ICC Enclosure")
    add_xls(to_xls, panel_name, ICC[building], SD[building])
    ICC_Tot += 1
  else:
    print("\trejected: " + str(ICC[building]))

  if HHW[building] != 0 and HHW[building] != "NONE":
    panel_name = (str(building) + " - HHW Enclosure")
    add_xls(to_xls, panel_name, HHW[building], SD[building])
    HHW_Tot += 1
  else:
    print("\trejected: " + str(HHW[building]))

print(SBC_Tot, TC_Tot, ICC_Tot, HHW_Tot)

df_panels = pd.DataFrame(to_xls)

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('panels.xlsx', engine='xlsxwriter')

# Convert the dataframe to an XlsxWriter Excel object.
df_panels.to_excel(writer, sheet_name='Sheet1', index=False)

# Close the Pandas Excel writer and output the Excel file.
writer.save()






