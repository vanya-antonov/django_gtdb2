from django.core.management.base import BaseCommand, CommandError

from gtdb2.models.user import User


class Command(BaseCommand):
    help = "Deletes entire database content"

    def handle(self, *args, **options):
        self.stdout.write(self.style.WARNING(
            "This will delete ALL(!) objects from the DB!"))
        answer = input("Please type 'DELETE' to confirm or press Enter to cancel: ")
        if answer != 'DELETE':
            self.stdout.write(self.style.ERROR(
                'Canceled'))
            return

        for user in User.objects.all():
            user.delete()
            self.stdout.write(self.style.SUCCESS(
                'Successfully deleted user "%s"' % user.name))
        self.stdout.write(self.style.SUCCESS('ALL DONE!'))

