Data files related to the phenotyping experiment in Spring/Summer 2024 measuring
rosette size, flowering time and seed size in three replicates each of 219 parental
and 427 F9 lines.

- flowering_time_raw_data
    Directory containing raw data files for flowering time.
    Files giving the date followed by '_ft10.csv' are from MementoDB collected 
    by Tom Ellis. They are successive backups of the same data, so are partially
    redundant.
    Other files are collected by Almudena Molla Morales, and give the plants that
    flowered on a single day. Three files contain 'Estimated', which are plants
    from replicate 3 that Almudena did not notice, because they were on a high
    shelf, and she estimated in retrospect.
    Note that I manually corrected entry 1137 9405x6100 rep2 3.38.H4 to 1137_3.38.H4_9405x6100 rep2_F8
- F9_seed_availability.csv
    A list of seed stocks processed by Tal Dahan in February 2024 to be measured
    for seed size. The seeds are the F9 progeny of the F8 plants that were
    genotyped.
- F10_seed_size_cataloge.csv
    A list of all the seed packets collected from F9 plants grown in the
    phenotyping experiment.
- format_ft10_from_AMM.R
    Almudena kept track of flowering while I was away, but she used a different
    app to log data. This file formats those data.
- randomise_ft10.R
    Script to generate a layout for the flowering time experiment.
    Note that I could not find seed for 8258, but I did find a seed package labelled
    9369x9434_rep2 that was not in the list, I substituted these by hand.
- randomise_expt.csv
    Output of randomise_ft10.R giving the full layout of the experiment in order
- randomise_expt.txt
    Output of randomise_ft10.R ordered by genotype, for printing labels.
- rosette_size_parents_F9s.csv
    Rosette size data for cohorts 2 and 3 measured with the app "Easy Leaf Area"
    Columns indicate tray, row and column positions, genotype, whether the plant
    is a cross or a parents, three size measurements from different angles, and
    the ID of the observer.
        TJE = Tom Ellis
        AMM = Almudena Molla Morales
        TDM = Tal Dahan-Meir
        RM = Rozi
        NS = Natasa