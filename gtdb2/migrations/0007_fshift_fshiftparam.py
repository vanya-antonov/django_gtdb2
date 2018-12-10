# Generated by Django 2.1.3 on 2018-12-06 17:48

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('gtdb2', '0006_auto_20181203_1835'),
    ]

    operations = [
        migrations.CreateModel(
            name='FShift',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('c_date', models.DateTimeField(auto_now_add=True)),
                ('name', models.CharField(max_length=255, unique=True)),
                ('descr', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('coord', models.IntegerField()),
                ('type', models.IntegerField()),
                ('origin', models.CharField(max_length=255)),
                ('start', models.IntegerField()),
                ('end', models.IntegerField()),
                ('strand', models.IntegerField()),
                ('seq', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.Seq')),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.User')),
            ],
            options={
                'db_table': 'fshifts',
            },
        ),
        migrations.CreateModel(
            name='FShiftParam',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('value', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('num', models.FloatField(blank=True, default=None, null=True)),
                ('data', models.TextField(blank=True, default=None, null=True)),
                ('parent', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.FShift')),
            ],
            options={
                'db_table': 'fshift_params',
            },
        ),
    ]
