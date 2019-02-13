from django.core.management.base import BaseCommand, CommandError

from gtdb2.models.user import User


class Command(BaseCommand):
    help = "Deletes entire database content"

    def handle(self, *args, **options):
        self.stdout.write("This will delete ALL(!) objects from the DB!")
        answer = input("Please type 'DELETE' to confirm or press Enter to cancel: ")
        if answer != 'DELETE':
            print('Canceled')
            return

        for user in User.objects.all():
            user.delete()
            self.stdout.write(self.style.SUCCESS(
                'Successfully deleted user "%s"' % user.name))

