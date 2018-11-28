from django.conf import settings
from django.core.management.base import CommandError

from gtdb2.models.user import User

from ._abstract import AbstractCommand


class Command(AbstractCommand):
    help = "create default user '%s'" % settings.GTDB_USER

    def handle(self, *args, **options):
        user_name = settings.GTDB_USER
        if User.objects.filter(name=user_name).count() > 0:
            raise CommandError("User '%s' already exists!" % user_name)
        user = User(name=user_name, descr='Default GeneTackDB user')
        user.save()
        self.print_success("Successfully created user '%s' with id '%s'" %
                           (user.name, user.id))

