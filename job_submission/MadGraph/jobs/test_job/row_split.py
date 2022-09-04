import pandas as pd
import math

def row_dropper(filepath, x, newfile):
	"""reads in a csv and drops x rows from the top"""
	old_df = pd.read_csv(filepath)
	print(old_df.tail(5))
	to_drop = list(range(0, x))
	new_df = old_df.drop(old_df.index[to_drop])
	print(new_df.head(5))

	new_df.to_csv(newfile, index=False)

#example of use of row_dropper
#row_dropper('/scratch/cb27g11/iridis_script/magellan01.csv', 11742, '/scratch/cb27g11/iridis_script/bq_mag_dropd.csv')

def csv_splitter(filepath, x, newnametag):
	"""reads in a csv located at 'filepath', splits this into 'x' csves with sizes as close to even as is possible. Names are assigned as 'newnametag'_00, 'newnametag'_01 and so on."""
	old_df = pd.read_csv(filepath)
	csv_len = float(len(old_df))
	x = float(x)
	print(csv_len/x)
	print(5.0/3.0)
	csv_size = int(math.ceil(csv_len/x))

	print("csv_size", csv_size)

	print("x", x)

	i = 1	

	while i !=x:
		print("Loop number: " + str(i))
		if i != x:
			
			print(" creating csv number: " + str(i))
			selected_rows = list(range(0, csv_size))
			print("Selected rows: " + str(selected_rows))
			
			new_df = old_df.iloc[selected_rows]
			print("new_df", new_df)
			
			print("Saving new_df...")
			new_df.to_csv(str(newnametag) + "_" + str(i) + ".csv", index=False)

			print("dropping saved rows from old_df. Before:  ")
			print(old_df)
			old_df = old_df.drop(old_df.index[selected_rows])
			print("old_df", old_df)

			i += 1

	if i == x:
		old_df.to_csv(str(newnametag) + "_" + str(i) + ".csv", index=False)


csv_splitter("/scratch/cb27g11/Looping_Bash/missing.csv", 2, "/scratch/cb27g11/THDM_T3PS_scanner/job_submission/MadGraph/jobs/test_job/split_csv_name")



