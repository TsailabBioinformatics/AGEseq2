'''
Helper class to store SNPs per Indel.  Structure:

{ INDEL: {
    SNP: freq,
    SNP: freq,
  },
}

Algorithm:
(Each line in the file simulates a coreEP list passed into the algorithm)

create master set of data
create set for indels per line
create list for snps per line

process line in file (we have to go through it all in order to backpair indel:snp)
    if the item is an indel via I or D being in the letter
        add to indels set to ensure uniqueness
        add to dictionary if it does not exist as a key
    else
        add to snps list

    if no indels -> return
    else
        for each snp in snps
            for each indel in indel
                if data[indel][snp] is not an entry -> set to 1
                if data[indel][snp] is an entry -> increment by 1

    empty indels set
    empty snps set




'''
from pathlib import Path
from pprint import pprint


class IndelSNPStore:

    def __init__(self):
        self.data = {}
        self.test_file_location = Path('./indels_snp_store_files/example_coreEP.txt')

        self.parse()

    def parse(self):
        with open(self.test_file_location, 'r') as f:
            # Processing the file to simulate data
            for line in f:
                items = line.strip().split(',')

                current_indels = set()
                current_snps = []
                for item in items:
                    item = item.strip()

                    if 'I' in item or 'D' in item:
                        current_indels.add(item)
                        # update dictionary to avoid second list
                        if item not in self.data.keys():
                            self.data[item] = dict()

                    else:
                        current_snps.append(item)

                if current_indels:
                    for snp in current_snps:
                        for indel in current_indels:
                            if snp not in self.data[indel].keys():
                                self.data[indel][snp] = 1
                            else:
                                self.data[indel][snp] += 1
        pprint(self.data)


newline = ""
if __name__ == '__main__':
    test = IndelSNPStore()
