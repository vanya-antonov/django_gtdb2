
import importlib
import logging
from pprint import pprint

from gtdb2.lib.baseutil import get_all_ids
from gtdb2.lib.command import AbstractCommand


class Command(AbstractCommand):
    help = "Create an object and call one of its methods."

    def add_arguments(self, parser):
        parser.add_argument('method_name', help="""method name, e.g.
                            'create_all_params' or '_make_param_taxonomy'""")
        all_cls = ['Org', 'Seq', 'Feat', 'Fshift', 'ChelataseOrg',
                   'ChelataseSeq', 'ChelataseFeat', 'ChelataseFshift']
        parser.add_argument('cls_name', metavar='cls_name', choices = all_cls,
                            help='valid values: ' + ', '.join(all_cls))
        parser.add_argument('input_ids', default='-',
                            help="""a single ID or a file name with a list of
                            IDs (no header) or read the IDs from STDIN
                            (default)""")
        parser.add_argument('method_args', nargs='*', default=[])

    def handle(self, *args, **options):
        super().handle(*args, **options)

        # Identify the source module
        if options['cls_name'].startswith('Chelatase'):
            module_name = 'chelatase_db.models'
        else:
            module_name = 'gtdb2.models'

        # Get the class object by its name from the module
        # https://stackoverflow.com/a/4821120/310453
        module = importlib.import_module(module_name)
        cls = getattr(module, options['cls_name'])

        # Create the object and execute the method
        for obj_id in get_all_ids(options['input_ids']):
            obj = cls.objects.get(pk=obj_id)
            method_args = options['method_args']
            logging.info(
                'Calling methond "%s" on object with id "%s" (%s) and args "%s"' %
                (options['method_name'], obj_id, obj, method_args))
            # Execute the method: https://stackoverflow.com/a/3521742/310453
            getattr(obj, options['method_name'])(*method_args)

