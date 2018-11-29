from django.core.management.base import BaseCommand


class AbstractCommand(BaseCommand):
    """This class adds common options and methods."""

    def print_success(txt):
        """A convenience method to pring messages in nice green color."""
        self.stdout.write(self.style.SUCCESS(txt))

