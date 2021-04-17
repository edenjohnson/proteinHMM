from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Create a mammals dictionary to write the common name of the species sequence in the data set.
mammals = {"NP_002015.1": "Human",
           "XP_003806864.1":"Bonobo",
           "XP_012354562.2":"Norther_white-cheeked_gibbon",
           "XP_024096084.1":"Sumatran_orangutan",
           "XP_005594814.1":"Crab-eating_macaque",
           "XP_011848408.1":"Drill",
           "XP_030789842.1":"Golden_snub-nosed_monkey",
           "XP_012322520.1":"Nancy_Ma's_night_monkey",
           "XP_008271672.1":"European_rabbit",
           "XP_004380984.1":"West_Indian_manatee",
           "XP_004777989.1":"Ferret",
           "XP_008570460.1":"Sunda_flying_lemur",
           "XP_004478534.1":"Nine-banded_armadillo",
           "XP_006869781.1":"Cape_golden_mole",
           "XP_004413934.1":"Pacific_walrus",
           "XP_032284155.1":"Harbor_seal",
           "XP_025720881.1":"Northern_fur_seal",
           "XP_036695669.1":"Blue_whale",
           "XP_001925767.1":"Wild_boar",
           "XP_022377002.1":"Sea_otter",
           "XP_039111949.1":"Striped_hyena",
           "XP_004000998.1":"Cat",
           "XP_036126371.1":"Velvety_free-tailed_bat",
           "XP_021538209.1":"Hawaiian_monk_seal",
           "XP_037367742.1":"Spanish_mole"}

# Create a birds dictionary to write the common name of the species sequence in the data set.
birds = {"XP_025935695.1":"Okarito_kiwi",
         "XP_009329388.1":"Adélie_penguin",
         "XP_009887558.1":"killdear",
         "XP_009273621.1":"Emperor_penguin",
         "XP_009932965.1":"Hoatzin",
         "XP_009580872.1":"Norther_fulmar",
         "XP_009899824.1":"Downy_woodpecker",
         "XP_010019670.1":"Kea",
         "XP_009875843.1":"Bar-tailed_trogon",
         "XP_009950999.1":"Cuckoo_roller",
         "XP_010079211.1":"Yellow-throated_sandgrouse",
         "NXF78372.1":"Tawny-throated_leaftosser",
         "NXM75263.1":"Silver-breasted_broadbill",
         "NWX03175.1":"Nicobar_pigeon",
         "NXF67291.1":"Black-and-white_Owl",
         "NXY67620.1":"Collared_pratincole",
         "NXO53151.1":"Limpkin",
         "NXN58551.1":"Black_skimmer",
         "NXJ49717.1":"Black_hawk-eagle",
         "NWV11778.1":"Satin_bowerbird",
         "NXB64047.1":"Apostlebird",
         "NWX76680.1":"Razorbill",
         "NXO90270.1":"Short-toed_treecreeper",
         "NXN95020.1":"Common_scimitarbill",
         "NXA17942.1":"Ibisbill"}

basic_training_file = "basic_training_set.fa"
basic_testing_file = "basic_testing_set.fa"
basic_set_fixed = "basic_set_fixed.fa"
related_training_file = "related_training_set.fa"
related_testing_file = "related_testing_set.fa"
related_set_fixed = "related_set_fixed.fa"

# relabeling the species to their common names for input into the tree
with open("aligned.fa", "r") as input_handle, open("aligned_fixed.fa", "w") as output_handle:
    sequences = SeqIO.parse(input_handle, "fasta")

    fixed_sequences = []

    for line in sequences:
        if line.id in mammals.keys():
            fixed_record = SeqRecord(line.seq, mammals[line.id])
            fixed_sequences.append(fixed_record)
        elif line.id in birds.keys():
            fixed_record = SeqRecord(line.seq, birds[line.id])
            fixed_sequences.append(fixed_record)

    SeqIO.write(fixed_sequences, output_handle, "fasta")

# Relabel the species in the basic data set and also create the training/testing subsets
with open("basic_aligned.fa", "r") as basic_input_handle:
    sequences = SeqIO.parse(basic_input_handle, "fasta")
    basic_training_data_set = []
    basic_testing_data_set = []
    basic_count = 0

    for line in sequences:
        record = SeqRecord(line.seq, mammals[line.id])
        if basic_count < 20: # 80% training
            basic_training_data_set.append(record)
        elif basic_count >= 20: # 20% testing
            basic_testing_data_set.append(record)
        basic_count += 1

    basic_data_set = basic_training_data_set + basic_testing_data_set

    SeqIO.write(basic_data_set, basic_set_fixed, "fasta")
    SeqIO.write(basic_training_data_set, basic_training_file, "fasta")
    SeqIO.write(basic_testing_data_set, basic_testing_file, "fasta")

# Relabel the species in the basic data set and also create the training/testing subsets
with open("related_aligned.fa", "r") as related_input_handle:
    sequences = SeqIO.parse(related_input_handle, "fasta")
    related_training_data_set = []
    related_testing_data_set = []
    related_count = 0

    for line in sequences:
        if related_count < 20: # 80% training
            record = SeqRecord(line.seq, birds[line.id])
            related_training_data_set.append(record)
        elif 20 <= related_count < 25: # 20% testing and do not include the Homo sapiens query sequence used in the alignment
            record = SeqRecord(line.seq, birds[line.id])
            related_testing_data_set.append(record)
        related_count += 1

    related_data_set = related_training_data_set + related_testing_data_set

    SeqIO.write(related_data_set, related_set_fixed, "fasta")
    SeqIO.write(related_training_data_set, related_training_file, "fasta")
    SeqIO.write(related_testing_data_set, related_testing_file, "fasta")

input_handle.close()
output_handle.close()
basic_input_handle.close()
related_input_handle.close()