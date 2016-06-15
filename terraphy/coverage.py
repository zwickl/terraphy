#!/usr/bin/env python
import sys
from random import sample, random

#should funcs here that use tkinterutils be moved there?
#something better needs to be done if this import fails
#from tkarg.tkinterutils import *


class CoverageMatrix(object):
    '''A class to hold information on the pattern of missing data in a matrix.
    '''
    def __init__(self, per_locus_taxon_sets=None):
        '''per_taxon_presence_absence is a dict holding locus presence absence, i.e. {taxon_label:[0, 1, 1]}
        per_locus_taxon_sets is a list of sets of taxa present for each locus
        '''
        self.taxa = set()
        self.per_taxon_presence_absence = {}
        self.per_locus_taxon_sets = []
        self.reference_taxon = None
        
        if per_locus_taxon_sets:
            self.fill_from_subsets(per_locus_taxon_sets)

    def fill_from_subsets(self, subsets):
        '''Each subset is a list of taxa present for a given locus
        '''
        self.per_locus_taxon_sets = [set(col) for col in subsets]
        for col in self.per_locus_taxon_sets:
            self.taxa |= col
        self.fill_taxa()

    def fill_random(self, num_taxa, num_loci, coverage, func, reference_taxon=False):
        '''Simulate a taxon coverage matrix by making the locus taxon sets
        '''

        taxa = [ 't%d' % t for t in xrange(num_taxa) ]
        self.taxa = set(taxa)
        subsets = [ [] for _ in xrange(num_loci) ]

        for x in xrange(num_loci):
            for y in xrange(num_taxa):
                #print func(x, y, coverage),
                if (reference_taxon and y == 0) or (random() < func(x, y, coverage)):
                    subsets[x].append(taxa[y])
        self.fill_from_subsets(subsets)
 
    def fill_random_locus_func(self, num_taxa, num_loci, func, min_coverage=0.0, max_coverage=1.0, reference_taxon=False):
        '''Simulate a taxon coverage matrix by making the locus taxon sets
        '''
        taxa = [ 't%d' % t for t in xrange(num_taxa) ]
        self.taxa = set(taxa)
        subsets = [ [] for _ in xrange(num_loci) ]

        for x in xrange(num_loci):
            coverage = -1.0
            while coverage < min_coverage or coverage > max_coverage:
                coverage = func()
            for y in xrange(num_taxa):
                #print func(x, y, coverage),
                if (reference_taxon and y == 0) or (random() < coverage):
                    subsets[x].append(taxa[y])
        self.fill_from_subsets(subsets)
 
    def fill_random_taxon_func(self, num_taxa, num_loci, func, min_coverage=0.0, max_coverage=1.0, reference_taxon=False):
        '''Simulate a taxon coverage matrix by making the locus taxon sets
        '''
        taxa = [ 't%d' % t for t in xrange(num_taxa) ]
        self.taxa = set(taxa)
        subsets = [ [] for _ in xrange(num_loci) ]
        
        cov_list = []
        for _ in taxa:
            while coverage < min_coverage or coverage > max_coverage:
                coverage = func()
            cov_list.append(coverage)
       
        print cov_list, sum(cov_list) / float(num_taxa)

        for x in xrange(num_loci):
            for y in xrange(num_taxa):
                if (reference_taxon and y == 0) or (random() < cov_list[y]):
                    subsets[x].append(taxa[y])
        self.fill_from_subsets(subsets)
    
    def __getitem__(self, index):
        '''Pull information out of the coverage matrix, indexing taxa by name or number'''
        if isinstance(index, str):
            if not self.per_taxon_presence_absence:
                self.fill_taxa()
            if index in self.per_taxon_presence_absence:
                return self.per_taxon_presence_absence[index]
            else:
                raise KeyError('Taxon %s not found in coverage matrix\n')
        else:
            assert(isinstance(index, int))
            if not self.per_locus_taxon_sets:
                self.fill_loci()
            if len(self.per_locus_taxon_sets) < index:
                return self.per_locus_taxon_sets[index]

    def fill_loci(self, recalculate=False):
        '''Fill the column attribute from the rows'''
        if self.per_locus_taxon_sets and not recalculate:
            sys.exit('columns already filled')
        elif not self.per_taxon_presence_absence:
            sys.stderr.write('WARNING: filling matrix columns from empty matrix rows')
        for sub_num in xrange(len(self.per_taxon_presence_absence[0])):
            self.per_locus_taxon_sets.append([tax[sub_num] for tax in self.per_taxon_presence_absence])

    def fill_taxa(self, recalculate=False):
        '''Fill the rows from the columns'''
        if self.per_taxon_presence_absence and not recalculate:
            sys.exit('taxa already filled')
        elif not self.per_locus_taxon_sets:
            sys.stderr.write('WARNING: filling taxa from empty columns')
        self.per_taxon_presence_absence = { tax:[] for tax in self.taxa }
        for col_num, col in enumerate(self.per_locus_taxon_sets):
            col_set = set(col)
            for tax in self.taxa:
                if tax in col_set:
                    self.per_taxon_presence_absence[tax].append(1)
                else:
                    self.per_taxon_presence_absence[tax].append(0)

    def calculate_statistics(self):
        if not self.per_taxon_presence_absence:
            self.fill_taxa()
        if not self.per_locus_taxon_sets:
            self.fill_loci()
        self.num_matrix_cells = len(self.per_taxon_presence_absence) * len(self.per_locus_taxon_sets)
        self.filled_matrix_cells = sum([sum(row) for row in self.per_taxon_presence_absence.values()])
        self.coverage_proportion =float(self.filled_matrix_cells) / self.num_matrix_cells
        #print self.num_matrix_cells, self.filled_matrix_cells, self.coverage_proportion

    def print_subset_vectors(self):
        for c in self.per_locus_taxon_sets:
            print ' '.join(c)

    def loci_per_taxon(self):
        if not self.per_taxon_presence_absence:
            self.fill_taxa()

        num_loci = len(self.per_locus_taxon_sets)

        counts = [0] * (num_loci + 1)
        for row in self.per_taxon_presence_absence.values():
            counts[sum(row)] += 1
        return counts
 
    def find_reference_taxon(self):
        '''Choose an arbitrary taxon that has full coverage across loci, if
        one exists.  Set the class member reference_taxon to the taxon label
        and return it'''
        if not self.per_taxon_presence_absence:
            self.fill_taxa()

        for taxon, coverage in self.per_taxon_presence_absence.iteritems():
            if sum(coverage) == len(self.per_locus_taxon_sets):
                self.reference_taxon = taxon
                break
        else:
            self.reference_taxon = None

        return self.reference_taxon

    def test_decisiveness(self):

        #count = True
        count = False

        tlist = list(self.taxa)

        if not self.find_reference_taxon():
            print 'NO reference taxon present (decisiveness test may not be accurate)' 
        else:
            print 'Reference taxon present' 
        
        num_tests = 0
        found = 0
        num_q_tests = 0
        #this is equivalent to what itertools.combinations would do, but for some
        #reason this is somewhat faster both under regular python and pypy
        for n1, t1 in enumerate(tlist, 1):
            for n2, t2 in enumerate(tlist[n1:], n1+1):
                for n3, t3 in enumerate(tlist[n2:], n2+1):
                    for n4, t4 in enumerate(tlist[n3:], n3+1):
                        q = {t1, t2, t3, t4}
                        num_q_tests += 1
                        #if num_q_tests % 10000 == 0:
                        #    print '%g' % num_q_tests
                        for col in self.per_locus_taxon_sets:
                            num_tests += 1
                            if q.issubset(col):
                                found += 1
                                break
                        else:
                            if count:
                                pass
                            else:
                                pass
                                #print 'FAILED',
                                return False
        #print  '\t'.join(['Tests: %g' % num_tests, 'Qs: %g' % num_q_tests, 'QsFound: %g' % found, 'Prop: %g' % (float(found)/num_q_tests)])
        return True
        
        num_tests = 0
        found = 0
        num_t_tests = 0
        if self.find_reference_taxon():
            tlist.remove(self.reference_taxon)
        for n1, t1 in enumerate(tlist, 1):
            for n2, t2 in enumerate(tlist[n1:], n1+1):
                for n3, t3 in enumerate(tlist[n2:], n2+1):
                    trip = {t1, t2, t3}
                    num_t_tests += 1
                    #if num_t_tests % 10000 == 0:
                    #    print '%g' % num_t_tests
                    for col in self.per_locus_taxon_sets:
                        num_tests += 1
                        if trip.issubset(col):
                            found += 1
                            break
                    else:
                        if count:
                            pass
                        else:
                            #print 'FAILED',
                            #pass
                            return False
        #print (True, '%g' % num_tests, '%g' % num_t_tests, '%g' % found, '%g' % (float(found)/num_t_tests))
        #print  '\t'.join(['Tests: %g' % num_tests, 'Ts: %g' % num_t_tests, 'TsFound: %g' % found, 'Prop: %g' % (float(found)/num_t_tests)])
        return True

    def draw_matrix_graphic(self, canvas, sorted_taxa, x_offset, y_offset, width, height):
        '''Draw a checkerboard representation of matrix coverage onto a tkinter canvas
        with the upper left corner being located at (x_offset, y_offset).
        The default is to choose the smaller of taxa/loci to be the rows
        So, this can transpose from the typical alignment view, which might be
        a little confusing.
        sorted_taxa can be used to reorder the rows, e.g. sorting by coverage density
        '''

        x_labels_height = 50
        y_labels_width = 50
        label_buffer = 5

        border_line_width = 2
        width -= border_line_width * 2 + y_labels_width + label_buffer
        height -= border_line_width * 2 + x_labels_height + label_buffer

        loci_on_x = True if (len(self.per_taxon_presence_absence) < len(self.per_locus_taxon_sets)) else False
        if not loci_on_x:
            num_x, num_y = len(self.per_taxon_presence_absence), len(self.per_locus_taxon_sets)
        else:
            num_y, num_x = len(self.per_taxon_presence_absence), len(self.per_locus_taxon_sets)

        x_box_size = max(1, width  / num_x)
        y_box_size = max(1, height / num_y)

        max_ratio = 8
        if x_box_size < y_box_size:
            y_box_size = min(y_box_size, x_box_size * max_ratio)
        else:
            x_box_size = min(x_box_size, y_box_size * max_ratio)

        #the amount of space actually used, due to rounding of box sizes
        actual_width = x_box_size * num_x + border_line_width * 2
        actual_height = y_box_size * num_y + border_line_width * 2
        
        x_rect_offset = x_offset + y_labels_width + label_buffer
        y_rect_offset = y_offset + x_labels_height + label_buffer
        x0, y0 = x_offset + border_line_width + x_labels_height + label_buffer, y_offset + border_line_width + y_labels_width + label_buffer
        x_loc, y_loc = x0, y0

        #rectangle around the checkerboard
        canvas.create_rectangle(x_rect_offset, y_rect_offset, x_rect_offset+actual_width+border_line_width, y_rect_offset+actual_height+border_line_width, width=2)

        #AXIS TITLES
        xtext, ytext = ('Locus', 'Taxon') if loci_on_x else ('Taxon', 'Locus')
        #using a Label for y axis to allow it to wrap in a single veritcal line, rotation is hard or not possible
        ylab = Label(canvas, text=ytext, wraplength=1)
        canvas.create_window((x_offset, y_rect_offset + actual_height / 2), window=ylab, height=ylab.winfo_reqheight(), width=ylab.winfo_reqwidth(), anchor='w')
        canvas.create_text(x_rect_offset + actual_width / 2, y_rect_offset - label_buffer - x_labels_height / 2, text=xtext, anchor='s') 

        #AXIS LABELS
        ideal_x_labels, ideal_y_labels = 10, 10
        if num_x < ideal_x_labels:
            x_labels = range(1, num_x+1)
        else:
            x_labels = range(1, num_x+1, num_x / ideal_x_labels)
        for lab in x_labels:
            canvas.create_text(x0 + x_box_size / 2 + x_box_size * (lab - 1), y_rect_offset - label_buffer, text=('%d' % lab), anchor='s') 
       
        if num_y < ideal_y_labels:
            y_labels = range(1, num_y+1)
        else:
            y_labels = range(1, num_y+1, num_y / ideal_y_labels)
        for lab in y_labels:
            canvas.create_text((x_rect_offset - label_buffer, y0 + y_box_size / 2 + y_box_size * (lab - 1)), text='%d' % lab, anchor='e') 

        #now the checkerboard
        if loci_on_x:
            #for tax, cov in self.per_taxon_presence_absence.items():
            for tax in sorted_taxa:
                cov = self.per_taxon_presence_absence[tax]
                for cell in cov:
                    if cell == 1:
                        canvas.create_rectangle(x_loc, y_loc, x_loc + x_box_size, y_loc + y_box_size, fill="blue", outline='blue')
                    else:
                        canvas.create_rectangle(x_loc, y_loc, x_loc + x_box_size, y_loc + y_box_size, fill='white', outline='white')
                    x_loc += x_box_size
                y_loc += y_box_size
                x_loc = x0
        else:
            for num, cov in enumerate(self.per_locus_taxon_sets, 1):
                for tax in sorted_taxa:
                    if tax in cov:
                        canvas.create_rectangle(x_loc, y_loc, x_loc + x_box_size, y_loc + y_box_size, fill="blue", outline='blue')
                    else:
                        canvas.create_rectangle(x_loc, y_loc, x_loc + x_box_size, y_loc + y_box_size, fill='white', outline='white')
                    x_loc += x_box_size
                y_loc += y_box_size
                x_loc = x0

    def draw_barplot(self, canvas, counts, x_offset, y_offset, width, height):
        '''Draw a barplot of something like loci-per-taxon. This actually just plots the passed
       in count list, and doesn't use any of the class data
        '''

        x_labels_height = 40
        y_labels_width = 40
        plot_buffer = 7

        max_x = len(counts) - 1
        max_y = max(counts)
        y_size_per_count = (height - x_labels_height - plot_buffer) / float(max_y)
        
        num_bars = len(counts)
        bar_spacing = 5
        bar_width = (width - plot_buffer - y_labels_width - (num_bars - 1) * bar_spacing) / num_bars

        x0 = x_offset + y_labels_width
        y0 = y_offset + height - x_labels_height

        canvas.create_line(x0, y_offset, x0, y0, width=2)
        canvas.create_line(x0, y0, x0 + width, y0, width=2)
        
        #AXIS LABELS
        #using a Label for y axis to allow it to wrap in a single veritcal line, rotation is hard or not possible
        ylab = Label(canvas, text='Frequency', wraplength=1)
        canvas.create_window((x_offset, height / 2), window=ylab, height=ylab.winfo_reqheight(), width=ylab.winfo_reqwidth(), anchor='w')
        canvas.create_text(x_offset + y_labels_width + (width - y_labels_width) / 2, y_offset + height - plot_buffer, text='Number of loci', anchor='c') 
        
        ideal_y_labels = 5
        if max_y < ideal_y_labels:
            y_labels = range(max_y+1)
        else:
            y_labels = range(0, max_y+1, max_y / ideal_y_labels)
        for lab in y_labels:
            canvas.create_text(x0 - plot_buffer, y0 - (lab * y_size_per_count), text=('%d' % lab), anchor='e') 

        print x0 - plot_buffer, y0 - (lab * y_size_per_count)

        ideal_x_labels = 10
        if max_x < ideal_x_labels:
            x_labels = range(max_x+1)
        else:
            x_labels = range(0, max_x+1, max_x / ideal_x_labels)

        x_loc = x0
        for num, count in enumerate(counts):
            if num in x_labels:
                canvas.create_text(x_loc + bar_width / 2, y0 + plot_buffer, text='%d' % num, anchor='n')
            if count:
                bar_height = count * y_size_per_count
                canvas.create_rectangle(x_loc, y0, x_loc + bar_width, y0 - bar_height, fill="blue", outline='blue')
            x_loc += bar_width + bar_spacing

    #def print_statistics(self, canvas, x_offset, y_offset, width_height):
    #    self.calculate_statistics()




