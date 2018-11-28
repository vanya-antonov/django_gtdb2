"""A collection of functions to handle manage.py commands."""

from gtdb2.lib.db import GeneTackDB


def add_org_gbk(gtdb, fn, args):
    for record in SeqIO.parse(fn, "genbank"):
        # Create new org, if needed
        org_id = Org.get_db_id_by_SeqRecord(gtdb, record)
        if org_id is None:
            org_id = Org.create_new_in_db_from_SeqRecord(gtdb, record, args.user_id)
            org = Org(gtdb, org_id)
            org.add_param('source_fn', fn)
            logging.info("New org '{}' has been created ({})".format(record.annotations['organism'], org_id))
       
        # Create new seq, if needed
        seq_id = Seq.get_db_id_by_ext_id(gtdb, record.id)
        if seq_id is None:
            org = Org(gtdb, org_id)
            seq_id = Seq.create_new_in_db_from_SeqRecord(gtdb, record, org, args.user_id)
            logging.info("New seq '{}' has been created ({})".format(record.id, seq_id))

class GeneTackDBManager:
    """A class to handle manage.py commands."""

    def __init__(self):
        self.db = GeneTackDB()

