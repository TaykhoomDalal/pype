from datetime import datetime

# url for searching ukbb data fields (if updated in future, supply new URL with all fields selected)
# I.E. select all fields in Stability, Strata, Item Type, Value Type, & leave 4 XXXX in srch= part of url
current_year = str(datetime.now().year)
DEFAULT_SHOWCASE_URL = 'https://biobank.ndph.ox.ac.uk/showcase/search.cgi?wot=0&srch=XXXX&sta0=on&sta1=on&sta2=on&sta3=on&sta4=on&str0=on&str3=on&str1=on&str2=on&fit0=on&fit10=on&fit20=on&fit30=on&fvt11=on&fvt21=on&fvt22=on&fvt31=on&fvt41=on&fvt51=on&fvt61=on&fvt101=on&yfirst=2000&ylast=' + current_year

# add more mappings here if you need to fix other covariates
dict_UKB_fields_to_names =  {'f.31.0.0': 'Sex', 
							'f.34.0.0': 'Year_of_birth', 
							'f.52.0.0': 'Month_of_birth',
							'f.53.0.0': 'Date_attended_center_0', 
							'f.53.1.0': 'Date_attended_center_1',
							'f.53.2.0': 'Date_attended_center_2', 
							'f.53.3.0': 'Date_attended_center_3',
							'f.54.0.0': 'Center_of_attendance_0',
							'f.54.1.0': 'Center_of_attendance_1',
							'f.54.2.0': 'Center_of_attendance_2',
							'f.54.3.0': 'Center_of_attendance_3',
							'f.21000.0.0': 'Ethnicity', 
							'f.21000.1.0': 'Ethnicity_1', 
							'f.21000.2.0': 'Ethnicity_2',
							'f.22001.0.0': 'Sex_genetic'}

# One hot encode ethnicity
dict_ethnicity_codes = {'1': 'Ethnicity.White', '1001': 'Ethnicity.British', '1002': 'Ethnicity.Irish',
						'1003': 'Ethnicity.White_Other',
						'2': 'Ethnicity.Mixed', '2001': 'Ethnicity.White_and_Black_Caribbean',
						'2002': 'Ethnicity.White_and_Black_African',
						'2003': 'Ethnicity.White_and_Asian', '2004': 'Ethnicity.Mixed_Other',
						'3': 'Ethnicity.Asian', '3001': 'Ethnicity.Indian', '3002': 'Ethnicity.Pakistani',
						'3003': 'Ethnicity.Bangladeshi', '3004': 'Ethnicity.Asian_Other',
						'4': 'Ethnicity.Black', '4001': 'Ethnicity.Caribbean', '4002': 'Ethnicity.African',
						'4003': 'Ethnicity.Black_Other',
						'5': 'Ethnicity.Chinese',
						'6': 'Ethnicity.Other_ethnicity',
						'-1': 'Ethnicity.Do_not_know',
						'-3': 'Ethnicity.Prefer_not_to_answer',
						'-5': 'Ethnicity.NA'}

# mapping of ethnicity category to smaller categories
dict_ethnicity_mapping = {'Ethnicity.White' : ['Ethnicity.White', 'Ethnicity.British', 'Ethnicity.Irish', 'Ethnicity.White_Other'],
						'Ethnicity.Mixed' : ['Ethnicity.Mixed','Ethnicity.White_and_Black_Caribbean','Ethnicity.White_and_Black_African', 'Ethnicity.White_and_Asian', 'Ethnicity.Mixed_Other'],
						'Ethnicity.Asian': ['Ethnicity.Asian', 'Ethnicity.Indian', 'Ethnicity.Pakistani', 'Ethnicity.Bangladeshi','Ethnicity.Asian_Other'],
						'Ethnicity.Black': ['Ethnicity.Black', 'Ethnicity.Caribbean', 'Ethnicity.African','Ethnicity.Black_Other'],
						'Ethnicity.Other': ['Ethnicity.Other_ethnicity', 'Ethnicity.Do_not_know','Ethnicity.Prefer_not_to_answer','Ethnicity.NA']}