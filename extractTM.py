def process(line, sequences):
    """ line with more than 20 residues are cut into 20 residues pieces

    :param line (string): a single line from input file, start with 1 or 0, followed by a sequence
    :param sequences (list): a list to add the processed lines

    :return: None
    """
    label = line[:2]
    seq = line[2:]
    num = len(seq)/20
    for i in range(num):
        sequences.append(label + seq[i*20:i*20+20])
    if len(seq) - num*20 > 2:
        sequences.append(label + seq[-20:])


if __name__ == '__main__':
    DATAFILE = 'YeastTM.dat'
    to_write = 'YeastTM20.dat'
    sequences = []
    with open(DATAFILE) as f:
        for line in f:
            if line != '':
                process(line.strip(), sequences)
    with open(to_write, 'w') as f:
        for seq in sequences:
            f.write('%s\n' % seq.encode('utf8'))