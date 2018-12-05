# Generated by Django 2.1.3 on 2018-11-28 17:29

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Org',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('c_date', models.DateTimeField(auto_now_add=True)),
                ('name', models.CharField(max_length=255)),
                ('descr', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('genus', models.CharField(max_length=255)),
                ('phylum', models.CharField(max_length=255)),
                ('kingdom', models.CharField(max_length=255)),
            ],
            options={
                'db_table': 'orgs',
            },
        ),
        migrations.CreateModel(
            name='OrgParam',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('value', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('num', models.FloatField(blank=True, default=None, null=True)),
                ('data', models.TextField(blank=True, default=None, null=True)),
            ],
            options={
                'db_table': 'org_params',
            },
        ),
        migrations.CreateModel(
            name='Seq',
            fields=[
                ('c_date', models.DateTimeField(auto_now_add=True)),
                ('name', models.CharField(max_length=255)),
                ('descr', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('id', models.CharField(max_length=255, primary_key=True, serialize=False)),
                ('type', models.CharField(max_length=255)),
                ('ext_id', models.CharField(max_length=255)),
                ('len', models.IntegerField()),
                ('org', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.Org')),
            ],
            options={
                'db_table': 'seqs',
            },
        ),
        migrations.CreateModel(
            name='SeqParam',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('value', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('num', models.FloatField(blank=True, default=None, null=True)),
                ('data', models.TextField(blank=True, default=None, null=True)),
                ('seq', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.Seq')),
            ],
            options={
                'db_table': 'seq_params',
            },
        ),
        migrations.CreateModel(
            name='User',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('c_date', models.DateTimeField(auto_now_add=True)),
                ('name', models.CharField(max_length=255)),
                ('descr', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('password', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('email', models.CharField(blank=True, default=None, max_length=255, null=True)),
            ],
            options={
                'db_table': 'users',
            },
        ),
        migrations.CreateModel(
            name='UserParam',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('value', models.CharField(blank=True, default=None, max_length=255, null=True)),
                ('num', models.FloatField(blank=True, default=None, null=True)),
                ('data', models.TextField(blank=True, default=None, null=True)),
                ('user', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.User')),
            ],
            options={
                'db_table': 'user_params',
            },
        ),
        migrations.AddField(
            model_name='seq',
            name='user',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.User'),
        ),
        migrations.AddField(
            model_name='orgparam',
            name='org',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.User'),
        ),
        migrations.AddField(
            model_name='org',
            name='user',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='gtdb2.User'),
        ),
    ]