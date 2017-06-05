from __future__ import division
import pysam
import os
import logging

from pipeline.tgpolyt_depths import get_reads_covering_region, get_tgpolyt_bases, count_tgpolyt, count_tg_t, hethom_call

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

test_dir = os.path.dirname(os.path.realpath(__file__))
bam_file = os.path.join(test_dir, "resources", "NA12878.bam")

'''
Test parts of the script tgpolyt_depths.py

'''

REGION = ['7', 117188661, 117188690,'TG-polyT']

def test_get_reads_covering_region():
    '''
    test get_reads_covering_region
    '''
    logger.info('testing get_reads_covering_region()')
    #count reads in specified region
    reads = get_reads_covering_region( bam_file, REGION)
    assert len(reads) == 21

    #count reads in region without reads
    reads = get_reads_covering_region( bam_file, ['8', 117188618, 117188756, 'TG-polyT'])
    assert len(reads) == 0


def test_count_tgpolyt():
    '''
    count TG polyT in
    :return:
    '''

    logger.info('testing count_tgpolyt()')

    tgpolyt_counts = count_tgpolyt(get_reads_covering_region( bam_file, REGION), REGION)

    assert ('(TG)11-7T', 16) in tgpolyt_counts[0].items()


def test_get_tgpolyt_bases():
    '''
    test get_tgpolyt_bases
    :return:
    '''
    logger.info('testing get_tgpolyt_bases()')

    bampysam = pysam.Samfile( bam_file, "rb")
    read_dict = {x.query_name: x for x in bampysam.fetch()}

    read = read_dict['HISEQ1:17:H947YADXX:1:2213:5880:62287']
    assert get_tgpolyt_bases(read,REGION[1], REGION[2]) == 'ATGTGTGTGTGTGTGTGTGTGTGTTTTTTTA'

    read = read_dict['HISEQ1:17:H947YADXX:2:1208:7897:88820']
    assert get_tgpolyt_bases(read, REGION[1], REGION[2]) == ''

    # sequence that possibly causes an issue when calling polyT
    read = read_dict['HISEQ1:17:H947YADXX:1:2104:17965:28539']
    assert get_tgpolyt_bases(read, REGION[1], REGION[2]) == 'ATGTGTGTGTGTTTGTGTGTGTGTTTTTTTA'


def test_count_tg_t_01():
    '''
    testing count_tg_t under different scenarios
    '''
    logger.info('testing count_tg_t()')
    # easy sequence
    variant = count_tg_t("ACTTTGTGTGTTTTTATGTG")
    assert variant == "(TG)3-5T"

    # ambigous sequence with (TG)5-2T and (TG)5-7T variants
    variant = count_tg_t("ATGTGTGTGTGTTTGTGTGTGTGTTTTTTTA")
    if variant != "(TG)5-7T":
        logger.warning('missing variant (TG)5-7T {}'.format(variant))
    assert variant == "(TG)5-2T"

    # polyT sequence ending in TG which should yield (TG)3-5T
    variant = count_tg_t("ACTTGTGTGTTTTTGTG")
    if variant != "(TG)3-5T":
        logger.warning('missing variant (TG)3-5T {}'.format(variant))
    assert variant == "(TG)3-4T"

    # precessing TG but not followed by polyT
    variant = count_tg_t("ACTTGAATGTGTTTTTGTG")
    if variant != "(TG)2-5T":
        logger.warning('missing variant (TG)2-5T {}'.format(variant))
    assert variant == "(TG)1-0T"


def test_hethom_call():
    '''
    testing different scenarios with hethom_call

    '''
    logger.info('testing hethom_call()')

    assert hethom_call(['11-7T'], [16 / 17]) == 'hom'

    assert hethom_call(['11-9T'], [16 / 17]) == 'hom alt'

    assert hethom_call(['11-9T'], [8 / 20]) == 'het'

    assert hethom_call(['11-9T'], [5 / 20]) == '--'

    assert hethom_call(['11-7T', '5-2T'], [16 /17, 1/17]) == 'hom/--'

    assert hethom_call(['11-9T', '5-2T'], [7 / 20, 8 / 20]) == 'het/het'

    assert hethom_call(['11-9T', '5-2T'], [4 / 20, 4 / 20]) == '--/--'

    assert hethom_call(['11-9T', '11-7T', '5-2T'], [10 / 30, 10 / 30, 10 / 30]) == 'het/het/het'

    assert hethom_call(['11-9T', '5-2T'], [16 / 20, 4 / 20], hom_freq_threshold = 0.75) == 'hom alt/--'

    assert hethom_call(['11-9T', '5-2T'], [7 / 20, 4 / 20], hom_freq_threshold =  0.35) == 'het/--'
