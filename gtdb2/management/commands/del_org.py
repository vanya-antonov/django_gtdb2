from django.core.management.base import BaseCommand, CommandError

from gtdb2.models.org import Org


class Command(BaseCommand):
    help = "Deletes org from database and all associated files from GTDB dir"

    def add_arguments(self, parser):
        parser.add_argument('id', type=int, help="database id")

    def handle(self, *args, **options):
        # https://docs.djangoproject.com/en/2.1/howto/custom-management-commands/
        org_id = options['id']
        try:
            Org.objects.get(id=org_id).delete()
        except Org.DoesNotExist:
            raise CommandError('Org "%s" does not exist' % org_id)

        self.stdout.write(self.style.SUCCESS('Successfully deleted org "%s"' %
                                             org_id))

