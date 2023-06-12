# drought_CV_eBird_reproducible

A code repository to accompany the manuscript "Drought influences habitat associations and abundances of birds in California's Central Valley."

### List of files the user needs to run the full workflow

-   'data/greatvalley_outline.shp' (Provided). This is a shapefile that delineates the Central Valley ecoregion

-   'data/fveg_WHRNHUM' - the Functional Vegetation LANDFIRE raster

-   'data/fveg_attributes.txt' - A .txt file describing the FVEG features

-   'data/CA_Aug2020_species_counts.csv' - One of two parts constituting eBird relational database, which contains information on species counts for all California eBird checklists as of August 2020. This csv has four columns:

    -   SAMPLING.EVENT.IDENTIFIER - a unique ID indicating the eBird checklist event to which this observation corresponds

    -   name_clean - The common name of the species, as given by eBird, with whitespace and punctuation removed

    -   SCIENTIFIC.NAME - The scientific name of the species, as given by eBird

    -   total_count - The count reported for the species (\>=1, or "X" indicating no count reported. 0s are implicit)

-   'data/CA_Aug2020_checklist_info.csv' - One of two parts constituting a relational database of eBird data, this file stores metadata for each checklist event. This csv has the following columns:

    -   SAMPLING.EVENT.IDENTIFIER - a unique ID indicating the eBird checklist event to which this observation corresponds. Used for joining to counts csv

    -   LATITUDE

    -   LONGITUDE	

    -   OBSERVATION.DATE

    -   TIME.OBSERVATIONS.STARTED	

    -   OBSERVER.ID - a unique ID grouping checklists submitted by the same eBird observer

    -   DURATION.MINUTES

    -   EFFORT.DISTANCE.KM

    -   PROTOCOL.CODE - the eBird sampling protocol category that the checklist followed

    -   EFFORT.AREA.HA

    -   ALL.SPECIES.REPORTED - a logical flag indicating whether the checklist is complete

    -   NUMBER.OBSERVERS

-   'data/ebird_protocol_codes.csv' - Table of eBird protocol definitions. Provided

-   'data/all_CA_breeding_birds.csv' - List of birds breeding in California that were considered for this analysis. Provided

-   'data/LandUseClassification_AcrossSites_FVEG.xlsx' - Table relating land cover features to habitat groups. Provided

-   'data/spei07.nc' - The SPEI raster file. Obtained from SPEIbase.

-   'intermediate/EVI_df.csv' - A data frame giving the average EVI in each grid cell in each season. See manuscript for details on processing.

-   'intermediate/NDWI_df.csv' - A data frame giving the average NDWI in each grid cell in each season. See manuscript for details on processing.

### Contact me

If you have any questions about the code in this repository, you can contact the corresponding author at ben.goldstein [at] berkeley.edu
