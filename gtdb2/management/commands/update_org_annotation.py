from django.core.management.base import BaseCommand, CommandError

from gtdb2.models.org import Org


class Command(BaseCommand):
    help = "Re-creates annotation for the org"

    def add_arguments(self, parser):
        parser.add_argument('id', type=int, help="database id")

    def handle(self, *args, **options):
        # https://docs.djangoproject.com/en/2.1/howto/custom-management-commands/
        org_id = options['id']
        try:
            # TODO: load() -- returns ChelataseOrg object if py_class prm is defined
            # org = Org.load(org_id)
            org = Org.objects.get(id=org_id)
        except Org.DoesNotExist:
            raise CommandError('Org "%s" does not exist' % org_id)

        # Remove any objects created by the previous method call
        anno_user = org.gtdb.get_or_create_annotation_user()
        ChelataseFeat.objects.filter(user=anno_user).delete()
        ChelataseFshift.objects.filter(user=anno_user).delete()

        org.create_annotation()

        self.stdout.write(self.style.SUCCESS(
            'Successfully updated annotation for org "%s"' % org_id))

