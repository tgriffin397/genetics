import gzip


def get_best():
    best = Activation([])
    for activation in Stat.activation_list:
        if not best.read_list:
            best = activation
        else:
            if activation.get_repeats() + activation.get_reads_with_n() < best.get_repeats() + best.get_reads_with_n():
                best = activation
    return best


class Stat:
    activation_list = []


class Activation:
    def __init__(self, read_list):
        self.read_list = read_list

    def get_reads_with_n(self):
        n = 0
        for read in self.read_list:
            n += read.reads_with_n
        return n

    def get_sequence_avg(self):
        avgs = 0
        for read in self.read_list:
            avgs += read.get_sequence_avg()
        return round(avgs / len(self.read_list))

    def get_total_reads(self):
        return len(self.read_list)

    def get_repeats(self):
        repeats = 0
        for read in self.read_list:
            repeats += read.total_repeats
        return repeats

    def get_n_percent(self):
        n_percentage = 0
        for read in self.read_list:
            n_percentage += read.get_n_percentage()
        return round(n_percentage / len(self.read_list), 2)

    def get_gc_percent(self):
        gc_percentage = 0
        for read in self.read_list:
            gc_percentage += read.get_gc_percent()
        return round(gc_percentage / len(self.read_list), 2)


class Read:
    def __init__(self, lines):
        self.read_list = []
        self.total_reads = 0
        self.total_repeats = 0
        self.reads_with_n = 0
        self.total_read_length = 0
        self.sequence_list = {}
        self.n_list = []
        self.lines = lines  # LIST
        self.sequence = self.lines[1]

        Read.log_sequence(self)  # Class method.

    def log_sequence(self):
        #Stat.read_list.append(self)
        #Stat.total_reads += 1
        self.total_read_length += len(self.sequence)
        # Checks for repeats and if any Ns are found in sequence.
        if 'N' in self.sequence:
            self.reads_with_n += 1
            self.n_list.append(self.sequence.count('N') / len(self.sequence))
        if self.sequence not in self.sequence_list:
            self.sequence_list.update({self.sequence: 1})
        else:
            self.sequence_list.update({self.sequence: self.sequence_list[self.sequence] + 1})
            self.total_repeats += 1

    def get_sequence_avg(self):
        val_length = len(self.sequence_list.values())
        key_length = 0
        for key in self.sequence_list.keys():
            key_length += len(key)
        return round(key_length / val_length)

    def get_gc_percent(self):
        counts_list = []
        total_counts = {
            'A': 0,
            'C': 0,
            'G': 0,
            'T': 0,
            'N': 0,
        }
        for char in self.sequence:
            if char in total_counts:
                total_counts.update({char: total_counts[char] + 1})
        gc_percent = round((total_counts['G'] + total_counts['C']) / sum(total_counts.values()) * 100, 2)
        counts_list.append(gc_percent)

        gc_val_mean = round(sum(counts_list) / len(counts_list), 2)
        return gc_val_mean

    def get_n_percentage(self):
        return round((sum(self.n_list) / len(self.n_list)) * 100, 2)


def process_read(file_contents):
    read_list = []
    # For every 4 lines, create a read represented by a list and add it to another list
    # representing total reads.
    count = 0
    for i in range(len(file_contents) // 4):
        read_list.append(Read((file_contents[count:count + 4])))
        count += 4

    return read_list


def results(activation):
    print(f'Reads in the file = {activation.get_total_reads()}:')
    print(f'Reads sequence average length = {activation.get_sequence_avg()}')
    print('')
    print(f'Repeats = {activation.get_repeats()}')
    print(f'Reads with Ns = {activation.get_reads_with_n()}')
    print('')
    print(f'GC content average = {activation.get_gc_percent()}%')
    print(f'Ns per read sequence = {activation.get_n_percent()}%')


def main_process():
    for _ in range(3):
        with gzip.open(input(), 'rt') as f:
            Stat.activation_list.append(Activation(process_read(f.read().split('\n'))))


main_process()
results(get_best())
