from pathlib import Path

def input(primer):
    var = []
    with open(primer, 'r') as inp:
        lines = inp.readlines()
        for line in lines:
            var.append(line.strip())
    return var

def reference(gene_position):
    id = []
    start_gene = []
    end_gene = []
    with open(gene_position, "r") as positions:
        lines = positions.readlines()
        for line in lines:
            column = line.strip().split(';')
            id.append(column[0])
            start_gene.append(int(column[1]))
            end_gene.append(int(column[2]))
    return id, start_gene, end_gene

def matches(id, start_gene, end_gene, fw, rv, fw_name, rv_name):
    count = 0

    genomes = []
    with open('/home/sagawa/Documentos/dnaA_project/biopy_blast/matches/genome_list', 'r') as input:
        lines = input.readlines()
        for line in lines:
            genomes.append(line.strip())

    dir_id = Path(f'/home/sagawa/Documentos/dnaA_project/biopy_blast/matches/output/id/{fw_name}x{rv_name}')
    if not dir_id.exists():
        dir_id.mkdir()

    with open(f'{dir_id}/{fw}x{rv}', 'w') as id_output:
        for sequence in genomes:
                if sequence.endswith('.fna'):
                        start = []
                        end = []
                        fw_table = f"/home/sagawa/Documentos/dnaA_project/biopy_blast/loop_dnaA/output/{fw}/{sequence}.xml"
                        rv_table = f"/home/sagawa/Documentos/dnaA_project/biopy_blast/loop_dnaA/output/{rv}/{sequence}.xml"
                        with open(fw_table, "r") as foward:
                            lines = foward.readlines()
                            for line in lines:
                                column = line.strip().split('	')
                                start.append(int(column[8]))
                                name = column[1]
                        with open(rv_table, "r") as reverse:
                            lines = reverse.readlines()
                            for line in lines:
                                column = line.strip().split('	')
                                end.append(int(column[9]))
                        neck = False
                        var = 0
                        for k in range(len(id)):
                            if name == id[k]:
                                var = k
                        for s in start:
                            if s > start_gene[var] and s < end_gene[var]:
                                if neck:
                                    break
                                for e in end:
                                    if s > e:
                                        if s - e < 25000:
                                            count = count + 1
                                            neck = True
                                            id_output.write(f'{id[var]}\n')
                                            break
                                    else:
                                        if e - s < 25000:
                                            count = count + 1
                                            neck = True
                                            id_output.write(f'{id[var]}\n')
                                            break
    return count

def main():

    fw = '400_region'
    rv = '1200_region'
    gene_position = '/home/sagawa/Documentos/dnaA_project/biopy_blast/matches/dnaA_position'

    input_fw = input(f'/home/sagawa/Documentos/dnaA_project/biopy_blast/matches/input/{fw}')
    input_rv = input(f'/home/sagawa/Documentos/dnaA_project/biopy_blast/matches/input/{rv}')
    id, start_gene, end_gene = reference(gene_position)

    count_fw = 0
    count_rv = 0

    with open(f'/home/sagawa/Documentos/dnaA_project/biopy_blast/matches/output/summary/{fw}x{rv}', 'w') as output: #Preparing the header
        output.write('###\t')
        for rv_header in input_rv:
            output.write(f'{rv_header}\t')

        for primer_fw in input_fw:
            output.write(f'\n{primer_fw}\t')
            count_fw +=1
            for primer_rv in input_rv:
                var = matches(id, start_gene, end_gene, primer_fw, primer_rv, fw, rv)
                output.write(f'{var}\t')
                count_rv += 1
                print(f'Major: {count_fw}/{len(input_fw)} Minor: {count_rv}/{len(input_rv)}')
            count_rv = 0

main()
# I. 200 x 400
# II. 200 x 1200
# III. 400 x 1200
