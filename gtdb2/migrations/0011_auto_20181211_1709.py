# Generated by Django 2.1.3 on 2018-12-11 17:09

from django.db import migrations, models
import django.utils.timezone


class Migration(migrations.Migration):

    dependencies = [
        ('gtdb2', '0010_auto_20181211_1704'),
    ]

    operations = [
        migrations.AlterField(
            model_name='fshift',
            name='c_date',
            field=models.DateTimeField(default=django.utils.timezone.now),
        ),
        migrations.AlterField(
            model_name='org',
            name='c_date',
            field=models.DateTimeField(default=django.utils.timezone.now),
        ),
        migrations.AlterField(
            model_name='seq',
            name='c_date',
            field=models.DateTimeField(default=django.utils.timezone.now),
        ),
    ]