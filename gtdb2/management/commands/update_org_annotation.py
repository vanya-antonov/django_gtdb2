
from django.core.management.base import BaseCommand, CommandError

from chelatase_db.models import ChelataseOrg
from gtdb2.lib.command import AbstractCommand
from gtdb2.models import Org, Feat, Fshift


class Command(AbstractCommand):
    help = "Re-creates annotation for the org"

    def add_arguments(self, parser):
        parser.add_argument('id', type=int, help="database id")

    def handle(self, *args, **options):
        super().handle(*args, **options)

        # https://docs.djangoproject.com/en/2.1/howto/custom-management-commands/
        org_id = options['id']
        try:
            # TODO: load() -- returns ChelataseOrg object if py_class prm is defined
            # org = Org.load(org_id)
            org = ChelataseOrg.objects.get(id=org_id)
        except Org.DoesNotExist:
            raise CommandError('Org "%s" does not exist' % org_id)

        # Remove old annotation
        anno_user = org.gtdb.get_or_create_annotation_user()
        Feat.objects.filter(user=anno_user, seq__org=org).delete()
        Fshift.objects.filter(user=anno_user, seq__org=org).delete()

        # Create new annotation from scratch
        org.create_annotation()
        org.create_all_params()

        self.stdout.write(self.style.SUCCESS(
            "Successfully updated annotation for org '%s' (id=%s)" %
            (org.name, org.id)))

