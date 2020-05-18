
import importlib
import re

from django.core.management.base import BaseCommand, CommandError

from gtdb2.lib.baseutil import get_all_ids
from gtdb2.models import Org, Feat, Fshift
from gtdb2.lib.command import AbstractCommand


class Command(AbstractCommand):
    help = "Export various sequences from the database."

    ALL_CLS_NAMES = ['Org', 'Seq', 'Feat', 'Fshift']

    def add_arguments(self, parser):
        parser.add_argument('prm_name', default=None,
                            help="""type of the sequence to export (e.g.
                            'seq_nt' or 'translation')""")
        all_cls_str = ', '.join(self.ALL_CLS_NAMES)
        parser.add_argument('cls_name', metavar='cls_name',
                            choices = self.ALL_CLS_NAMES,
                            help='valid values: ' + all_cls_str)
        parser.add_argument('input_ids', nargs='?', default='-',
                            help="""a single ID or a file name with a list of
                            IDs (no header) or read the IDs from STDIN
                            (default)""")

        parser.add_argument(
            '--seqname', metavar='FMT', default='__INPUT_ID__',
            help="""[__INPUT_ID__] seqname format.
            <FMT> may include the following special names: __INPUT_ID__,
            __ORG_ID__, __ORG_NAME__, __ORG_KINGDOM__, __ORG_GENUS__,
            __ORG_PHYLUM__, __ORG_DIR_NAME__,
            __SEQ_ID__, __SEQ_NAME__, __SEQ_DESCR__,
            __FEAT_ID__, __FEAT_NAME__, __FEAT_DESCR__, __FEAT_PRM_GENE__,
            __FSHIFT_ID__, __FSHIFT_NAME__, __FSHIFT_DESCR__, etc""")
        parser.add_argument('--upper', action='store_true',
                            help='capitalize all sequence letters')
        parser.add_argument('--add_frame', action='store_true',
                            help='indicate the frame by _')

    def handle(self, *args, **options):
        super().handle(*args, **options)

        # Get the class object by its name from the module
        # https://stackoverflow.com/a/4821120/310453
        module_name = 'gtdb2.models'
        module = importlib.import_module(module_name)
        cls = getattr(module, options['cls_name'])

        all_ids = get_all_ids(options['input_ids'])
        prm_name = options['prm_name']
        for db_id in all_ids:
            obj = cls.objects.filter(pk=db_id).first()

            if obj is None:
                self.stderr.write("No object %s with id=%s" % (cls, db_id))
                continue

            if prm_name not in obj.prm:
                self.stderr.write(
                    "Object %s (%s, id=%s) doesn't have prm '%s'" %
                    (obj, cls, db_id, prm_name))
                continue

            seq = obj.prm[prm_name]
            if options['upper']:
                seq = seq.upper()
            if options['add_frame']:
                seq = re.compile(r'(...)').sub(r'\1_', seq).strip('_')

            seqname = _get_seqname_for_fmt(options['seqname'], obj)

            # About using self.stdout.write() instead of print()
            # https://docs.djangoproject.com/en/2.2/howto/custom-management-commands/
            self.stdout.write(">%s\n%s" % (seqname, seq))


def _get_seqname_for_fmt(fmt, obj):
    """ 
    '__INPUT_ID__|__SEQ_ID__|__ORG_NAME__|__FEAT_NAME__'
    =>
    '3514|NZ_GG770557.1|Mycobacterium parascrofulaceum|HMPREF0591_RS22050'
    """
    fmt = re.compile('__INPUT_ID__').sub(str(obj.id), fmt)

    if type(obj) == Org:
        fmt = _fmt_process_org(fmt, obj)
    elif type(obj) == Feat:
        fmt = _fmt_process_feat(fmt, obj)
        fmt = _fmt_process_seq(fmt, obj.seq)
        fmt = _fmt_process_org(fmt, obj.seq.org)
    elif type(obj) == Fshift:
        fmt = _fmt_process_fshift(fmt, obj)
        fmt = _fmt_process_seq(fmt, obj.seq)
        fmt = _fmt_process_org(fmt, obj.seq.org)
    else:
        raise TypeError("Unknown type = '%s'", type(obj))

    return fmt

def _fmt_process_feat(fmt, info):
    fmt = re.compile('__FEAT_ID__').sub(str(info.id), fmt)
    fmt = re.compile('__FEAT_NAME__').sub(str(info.name), fmt)
    fmt = re.compile('__FEAT_DESCR__').sub(str(info.descr), fmt)
    fmt = re.compile('__FEAT_TYPE__').sub(str(info.type), fmt)
    fmt = re.compile('__FEAT_PRM_GENE__').sub(str(info.prm.get('gene', '')), fmt)
    return fmt

def _fmt_process_fshift(fmt, info):
    fmt = re.compile('__FSHIFT_ID__').sub(str(info.id), fmt)
    fmt = re.compile('__FSHIFT_NAME__').sub(str(info.name), fmt)
    fmt = re.compile('__FSHIFT_DESCR__').sub(str(info.descr), fmt)
    fmt = re.compile('__FSHIFT_COORD__').sub(str(info.coord), fmt)
    fmt = re.compile('__FSHIFT_LEN__').sub(str(info.len), fmt)
    return fmt

def _fmt_process_seq(fmt, info):
    fmt = re.compile('__SEQ_ID__').sub(str(info.id), fmt)
    fmt = re.compile('__SEQ_NAME__').sub(str(info.name), fmt)
    fmt = re.compile('__SEQ_DESCR__').sub(str(info.descr), fmt)
    return fmt

def _fmt_process_org(fmt, info):
    # use str(...) everywhere to avoid errors when attribute is None
    fmt = re.compile('__ORG_ID__').sub(str(info.id), fmt)
    fmt = re.compile('__ORG_NAME__').sub(str(info.name), fmt)
    fmt = re.compile('__ORG_GENUS__').sub(str(info.genus), fmt)
    fmt = re.compile('__ORG_PHYLUM__').sub(str(info.phylum), fmt)
    fmt = re.compile('__ORG_KINGDOM__').sub(str(info.kingdom), fmt)
    fmt = re.compile('__ORG_DIR_NAME__').sub(str(info.prm.dir_name), fmt)
    return fmt

