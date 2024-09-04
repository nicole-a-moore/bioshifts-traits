The BioShifts.zip provides the full geodatabase 'BioShifts: a global geodatabase of climate-induced species redistribution over land and sea'.	
"It includes (1) the species range shifts extracted from the literature (BioShifts.csv), (2) associated reference list (References.csv), and (3) an ArcGIS File Geodatabase (Study_Areas.gdb) containing one shapefile for each of the study area (Source)."	
	
We kindly ask the users to cite Lenoir et al. as well as the BioShifts dataset in any published material produced using these data:	
"Lenoir J., Bertrand, Comte L., R., Murienne J., Bourgeaud L., Hattab T., Grenouillet G.  (2020) Species better track climate warming in the oceans than on land. TO BE COMPLETED"	
"Comte L., Grenouillet G., Bertrand, R., Murienne J., Bourgeaud L., Hattab T., Lenoir J. (2020) BioShifts: a global geodatabase of climate-induced species redistribution over land and sea. Figshare, 10.6084/m9.figshare.7413365."	
	
BioShifts.csv: species range shifts database where each row corresponds to a species range shift estimate extracted from the litterature	
Kingdom	Taxonomic kingdom of the species for which the range shift has been estimated (according to NCBI; September 2019) 
Phylum	Taxonomic phylum of the species for which the range shift has been estimated (according to NCBI; September 2019)
Order	Taxonomic order of the species for which the range shift has been estimated (according to NCBI; September 2019)
Class	Taxonomic class of the species for which the range shift has been estimated (according to NCBI; September 2019)
Family	Taxonomic family of the species for which the range shift has been estimated (according to NCBI; September 2019)
Genus	Taxonomic genus of the species for which the range shift has been estimated (according to NCBI; September 2019)
Species	Validated species name of the species for which the range shift has been estimated (according to NCBI; September 2019)
Hemisphere	Hemisphere of the study: North or South
Ecosystem	Ecosystem of the study: Marine or terrestrial (note that freshwater species have been classified as terrestrial)
Gradient	Spatial gradient along which the species range shift has been estimated: Elevation or Latitude
Position	"Range shift parameter for which the range shift has been estimated: Leading edge, Trailing edge, Centroid"
ShiftR	Range shift estimate standardized by the length of the study period over which the range shift has been estimated
Unit	Unit of the range shift estimate: m/year or km/year
PrAb	Type of the data used to estimate the range shift: Occurrence (OCCUR) or Abundance (ABUND)
Sampling	"Characteristics of the time periods over which the range shift has been estimated: Continuous (CONT; yearly data or less), Irregular (IRR; multiple periods irregularly distributed), Multiple (MULT; multiple periods regularly distributed), Two (TWO; two periods)"
Grain	"Spatial grain of the data used to estimate the range shift: fine (FINE; spatial resolution lower than 10 km), coarse (COARSE; spatial resolution greater than 100 km), or medium (MEDIUM; intermediate situations)"
Signif	Whether the significance of the range shift estimate has been tested in the original publication (note that we do not report the result of such test in the database): YES or NO
Quality	"The quality of the approach used to estimate the range shift: LOW (no pre-processing of the raw data), BALANCED (data cleaning or resampling procedures were carried out to quantify the range shift on a balanced dataset), MODELED (range shifts were quantified species ranges modelled using  species distribution models calibrated independently for each time period); RESURVEYED (range shifts were quantified from paired designs such as permanent plots)"
Start	The first year of the temporal period over which the range shift has been estimated
Reference	The reference to the original publication
Source	Unique identifier of the spatial polygons delinating the study areas (provided in the geodatabase)
