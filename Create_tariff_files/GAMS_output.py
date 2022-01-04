import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime as dt
from datetime import timedelta as td
import sys 

import interpolate_matrix_test as interpolate_helpers


###############################
	# 1. Set parameters
	# 2. READ IN IN FILES
	# 3. SET DATA OBJECTS
	# 4. WRITE OUT RESULTS TO REQUIRED TXT AND CSV FILES REQUIRED FOR OPTIMISATION
###############################

###############################
	#TODO:
		# 1. Change .xlsx files to .csv files
		# 2. Fix folder reference files 
		# 3. interpolate_helpers.output_write_df writes multiple dfs when > 1 columns in 'input_df'
###############################




#### SET PARAMETERS 
year_length = 8760
interval_length = 1
repeat_ren_signals = 20
repeat_other_signals = 20

if interval_length == 1:
	interval_len = "_hourly"
elif interval_length == 4:
	interval_len = "_15min"
elif interval_length == 12:
	interval_len = "_5min"

month_vec = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

dir2 = 'C:/Users/w47147/misc_code/RODeO-master/RODeO-master/Create_tariff_files/Data_files/'

dir1 = dir2 + "CSV_data/"
dir0 = dir2 + "TXT_files/" 




#### READ IN IN FILES
num1int = pd.read_csv(dir2 + "CSV_data/e_tou_8760.csv")
num1A = pd.read_csv(dir2 + "CSV_data/e_prices.csv")

num6A = pd.read_csv(dir1 + "d_tou_prices.csv")

numlBint = pd.read_excel(dir1 + "GAMS_Energy_Sale.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num1B = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_Energy_Purchase.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num1BB = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_AS.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num2 = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_FUEL.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num3 = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_renewables.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num4 = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "linear")

numlBint = pd.read_excel(dir1 + "GAMS_additional_load_profiles.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num9 = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_product_price.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num10A = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_product_consumed.xlsx", skiprows = 1, sheet_name = 'Sheet1')
for column in  numlBint.columns: 
	if column not in ["Date", "Interval"]:
		numlBint[column] = numlBint[column]/sum(numlBint[column]) * year_length/24
num11A = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")


numlBint = pd.read_excel(dir1 + "GAMS_Max_input_cap.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num12A = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")

numlBint = pd.read_excel(dir1 + "GAMS_Max_output_cap.xlsx", skiprows = 1, sheet_name = 'Sheet1')
num13A = interpolate_helpers.interpolate_matrix(numlBint, year_length, interval_length, "repeat")


num1int_ = num1int.copy()
for column in num1int_.columns:
	for price_integer in num1A.index:
		num1int_.loc[num1int_[column] == price_integer, column] = num1A.loc[price_integer][str(float(column))]

num1 = interpolate_helpers.interpolate_matrix(num1int, year_length, interval_length, "repeat")
num1_ = interpolate_helpers.interpolate_matrix(num1int_, year_length, interval_length, "repeat")
num1A1 = num1_.fillna(0).copy()
num5 = pd.read_csv(dir2 + "CSV_data/d_flat_prices.csv")

num8 = pd.read_csv(dir1 + "fixed_charge.csv")

numA =  pd.read_csv(dir1 + "tariff_property_list.csv")
Scenarios1 = numA['label']

num5int = pd.read_csv(dir2 + "CSV_data/d_flat_prices.csv")
num6int = pd.read_csv(dir2 + "CSV_data/d_tou_8760.csv")
num6int2 = num6int.copy()
num6Aint = pd.read_csv(dir2 + "CSV_data/d_tou_prices.csv")

#### Update TOU pricing
# lookup rates for TOU charges
num6int2 = num6int.copy()
# Convert integer values to price values from URDB
# TODO: clean this mess up
for column in num6int.columns:
	for price_integer in num6Aint.index:
		num6int2.loc[num6int2[column] == price_integer, column] = num6Aint.loc[price_integer][str(float(column))]

num6 = interpolate_helpers.interpolate_matrix(num6int2, year_length, interval_length, "repeat")
num6B = interpolate_helpers.interpolate_matrix(num6int, year_length, interval_length, "repeat")







#### SET DATA OBJECTS
filenames = [x.replace(" ", "_") + ".txt" for x in Scenarios1]
file_names = [x.replace(" ", "").replace("%", "") + ".txt" for x in filenames]

# Create the starting date as a `datetime` object.
start = dt(1900, 1, 1, 0, 0, 0)
# List initialiser.
result = [start]

# Build a list of datetime objects for each hour of the year.
for i in range(1, 8760):
    start += td(seconds=3600)
    result.append(start)

# Initialise a DataFrame data structure.
df = pd.DataFrame({'dates': result})
# Add each column by extracting the object of interest from the datetime.
df['8760'] = df.index+1
df['month'] = df['dates'].dt.month

GAMS_string_df = df.groupby("month").agg({"8760":["min", "max"]})
GAMS_string_df["GAMS_string"] = "/"+ GAMS_string_df["8760"]["min"].astype(str) + "*" + GAMS_string_df["8760"]['max'].astype(str) +"/"
GAMS_string_df["Month_to_hour"] = GAMS_string_df.index.astype(str) + "." + GAMS_string_df["8760"]["min"].astype(str) + "*" + GAMS_string_df["8760"]['max'].astype(str)
month_to_hour_string = "/ " + ",".join(GAMS_string_df["Month_to_hour"]) + " /"




######################
# START RE-WRITE -> SHOULD BE ABLE TO AVOID FOR LOOPS AND LEVERAGE DATA STRUCTURE BETTER
######################
# Format month string for GAMS optimisation
utility_month_tranch = { }
for column in num6B.columns:
	utility_month_tranch[column] = {}
	
	for month in month_vec:# GAMS_string_df.index:
		utility_month_tranch[column][month] = {}
		
		for tranch in num6B[column].unique():
			utility_month_tranch[column][month][tranch] = [ ]
	
	for index_,row in GAMS_string_df.iterrows():
		for i in range(row['8760']['min'], row['8760']['max']):
			if (num6B[column][i] != num6B[column][i-1]):
				temp_string = temp_string + "*" + str(i)
				utility_month_tranch[column][month_vec[int(index_)-1]][num6B[column][i-1]].append(temp_string)
				temp_string = str(i+1)
			if (i ==row['8760']['min']):
				temp_string = str(row['8760']['min'])
				
			if (i ==row['8760']['max']-1):
				temp_string = temp_string + "*" + str(i+1)
				utility_month_tranch[column][month_vec[int(index_)-1]][num6B[column][i-1]].append(temp_string)

for column in num6B.columns:
	for month in month_vec:# GAMS_string_df.index:
		for tranch in num6B[column].unique():
			utility_month_tranch[column][month][tranch] = "/" + ",".join(utility_month_tranch[column][month][tranch]) + "/"

# Create TOU charge to be written to GAMS optimisation file
string_list = [ ]
for column in num1int.columns:
	num1int['e_tou_lead' + column] = num1int[column].shift(1)
	num1int_filtered = num1int[num1int[column] != num1int['e_tou_lead' + column]]
	
	num1int_filtered['hour'] = num1int_filtered.index + 1
	num1int_filtered['hour_lagged'] = num1int_filtered['hour'].shift(-1, fill_value=0) - 1
	num1int_filtered[column] += 1
	
	num1int_filtered['tou_energy_charge_str'] = num1int_filtered[column].astype(str) + "." + num1int_filtered['hour'].astype(str) + "*" + num1int_filtered['hour_lagged'].astype(int).astype(str) 
	temp_string = num1int_filtered['tou_energy_charge_str'].to_list()
	
	temp_string[-1] = temp_string[-1].replace("nan","8760")
	
	string_list.append(temp_string)
print(num1int_filtered)
######################
# END RE-WRITE -> SHOULD BE ABLE TO AVOID FOR LOOPS AND LEVERAGE DATA STRUCTURE BETTER
######################
	
	
	
	


#### WRITE OUT RESULTS TO REQUIRED TXT AND CSV FILES
Inputs1 = ['elec_purchase_price(interval)','elec_sale_price(interval)','nonspinres_price(interval)','regdn_price(interval)','regup_price(interval)','spinres_price(interval)','meter_mnth_chg(interval)','Fixed_dem(months)','Timed_dem(timed_dem_period)','TOU_energy_prices(TOU_energy_period)']

tariff_files = "Tariff_files"
Path(tariff_files).mkdir(parents=True, exist_ok=True)
gams_string = GAMS_string_df["GAMS_string"].to_list()

# loop over available tariffs and write out to tariff.txt and GAMS optimisation files
for index, file_name in enumerate(filenames):
	txt_write_list = [ ]
	txt_write_list.append("$oneempty")
	txt_write_list.append("\nSet \n")
	
	# Add Month intervals to the string list
	# e.g. 	Month_Jan(interval)	hours	/1*744/
	for month_gamsstring in zip(month_vec, gams_string):
		txt_write_list.append(f'\tMonth_{month_gamsstring[0]}(interval)\thours\t{month_gamsstring[1]}\n')
	
	txt_write_list.append("\n")
	
	# TODO: add multiple tariff file lookup if applicable, check against MATLAB code
	# Add Month_tiers to the string list
	# e.g. 	Jan_4(interval)	hours	/33*46
	for month in month_vec:
		for tariff_num in sorted(utility_month_tranch[str(index+1)][month].keys()):
			txt_write_list.append(f'\t{month}_{tariff_num+1}(interval)\thours\t{utility_month_tranch[str(index+1)][month][tariff_num]}\n')

	# Add month to interval relation
	# e.g. 	month_interval(months,interval)	month to interval relation	/ 1.1*744,
	txt_write_list.append(f'\n\tmonth_interval(months,interval)\tmonth to interval relation\t{month_to_hour_string}\n \n')
	
	# Add electricity time of use bins
	# e.g. 	elec_TOU_bins(TOU_energy_period,interval)	Energy TOU bins for retail electricity / 5.1*32
	txt_write_list.append(f'\telec_TOU_bins(TOU_energy_period,interval)\tEnergy TOU bins for retail electricity / {",".join(string_list[index])} /')
	txt_write_list.append("\n")

	txt_write_list.append(f'\nparameter	{Inputs1[0]}	/ \n')	
	# Janky code; general clean
	num1A1['hour'] = num1A1.index + 1
	num1A1["1"] = num1A1["1"] * 1000
	num1A1['string'] = num1A1['hour'].astype(str) + "\t" + num1A1['1'].round(1).astype(str)
	string = "\n".join(num1A1['string'].to_list())
	string = string + "/;"
	txt_write_list.append(string)
	

	# iterate through five columns of Inputs and write out zero data frame. 
	for input in Inputs1[1:6]:
		txt_write_list.append(f'\nparameter	{input}	/ \n')
		
		#TODO: follow up with Josh on rationale behind this not sure what these zero dataframes are created; are data structures used for GAMS write/decision variables? 
		temp_zero_df = num1A1.copy()
		temp_zero_df['hour'] = temp_zero_df.index + 1
		temp_zero_df['string'] = temp_zero_df['hour'].astype(str) + "\t" + str(0)
		string = "\n".join(temp_zero_df['string'].to_list())
		string = string + "/;"
		txt_write_list.append(string)

	txt_write_list = interpolate_helpers.append_tariff(num8, txt_write_list, Inputs1[6])
	txt_write_list = interpolate_helpers.append_tariff(num5int.astype(int), txt_write_list, Inputs1[7])
	num6A = num6A * 1000
	txt_write_list = interpolate_helpers.append_tariff(num6A.round(0).astype(int), txt_write_list, Inputs1[8])
	num1A = num1A * 1000
	txt_write_list = interpolate_helpers.append_tariff(num1A.round(1), txt_write_list, Inputs1[9])
	
	with open(tariff_files + "/test_" + file_name, 'w') as f:
		f.writelines(txt_write_list)
	f.close()


	# Ancillary_services file
	num2.drop("Date", axis = 1, inplace = True)
	num2.columns = ['Interval','Nonspinning reserve','Regulation Down','Regulation Up','Spinning Reserve']
	num2.to_csv(dir0 + "Ancillary_services_" + interval_len + "_test.csv", index = False)
	
	# Simplify code into a function from 323 -> 386 to avoid repeating code 
	# To follow up; output for Energy Sales Price will have column names of integer ranking
	interpolate_helpers.output_write_df(num1B, dir0 + "Energy_sale_prices_" + column + interval_len + "_test_.csv")
	interpolate_helpers.output_write_df(num1BB, dir0 + "Energy_purchase_prices_" + column + interval_len + "_test_.csv")
	interpolate_helpers.output_write_df(num3, dir0 + "NG_" + column + interval_len + "_test_.csv")

	
	for column in num4:
		if column == "Interval":
			continue
		num4[column] = num4[column]/max(num4[column])	
	interpolate_helpers.output_write_df(num4, dir0 + "renewable_profiles_" + column + interval_len + "_test_.csv", repeat_ren_signals, 0)	
	interpolate_helpers.output_write_df(num9, dir0 + "Additional_load_" + column + interval_len + "_test_.csv", repeat_ren_signals, 0)	
	interpolate_helpers.output_write_df(num10A, dir0 + "Product_price_" + column + interval_len + "_test_.csv", repeat_ren_signals, 1)
	interpolate_helpers.output_write_df(num11A, dir0 + "Product_consumption_" + column + interval_len + "_test_.csv", repeat_ren_signals, 0)
	interpolate_helpers.output_write_df(num12A, dir0 + "Max_input_cap_" + column + interval_len + "_test_.csv", repeat_ren_signals, 1)
	interpolate_helpers.output_write_df(num13A, dir0 + "Max_output_cap_" + column + interval_len + "_test_.csv", repeat_ren_signals, 1)
	
	temp_df = pd.DataFrame({'current_interval':[-1],'next_interval':[1],'current_storage_lvl':[0.5],'current_monthly_max':[0.8],'max_interval':[""]})
	temp_df.to_csv(dir0 + "controller_input_values_test.csv", index = False)
	
	temp_df = pd.DataFrame()
	temp_df["Year"] = range(1,21)
	temp_df["MARCS_depreciation_percentage"] = [0.2,0.32,0.192,0.115,0.115,0.058,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
	temp_df.to_csv(dir0 + "MACRS_depreciation_schedule_test.csv", index = False)
	
	
	temp_df = pd.DataFrame()
	temp_df["Year"] = range(1,13)
	temp_df["MARCS_depreciation_percentage"] = [27.03,27.45,27.93,27.77,27.19,27.1,27.15,27.46,27.39,28.19,27.93,28.44]
	temp_df.to_csv(dir0 + "NSCR_test.csv", index = False)
