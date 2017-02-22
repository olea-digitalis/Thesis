def remove_columns(infile,outfile,k=2):
    with open(infile,'r') as in_f:
        first_line = in_f.readline()
        
        with open(outfile,'w') as out_f:
            out_f.write(first_line)
            line_count = 0
            none_count = 0
            for line in in_f:
                s = line
                cols = s.split('\t')
                line_count += 1
                if cols[0] == 'None':
                    none_count += 1
                    continue
                """
                    if cols[2] != 'None':
                        cols[0] = cols[2]
                    else:
                        continue
                """
                out_s = ""
                for c in range(k):
                    out_s += cols[c] + '\t'
                out_s = out_s.strip()
                out_s += '\n'
                out_f.write(out_s.replace(';','|'))
    print("Line Count:",line_count)
    print("None Count:",none_count)

remove_columns("bio_data_big_frag.txt","big_frag_fixed_delims.txt")
