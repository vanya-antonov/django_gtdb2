# Generated by Django 2.1.3 on 2019-02-04 17:49

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('gtdb2', '0011_auto_20190116_1751'),
    ]

    operations = [
        migrations.AlterField(
            model_name='feat',
            name='origin',
            field=models.CharField(choices=[('annotation', 'GenBank annotation'), ('genetack', 'GeneTack prediction'), ('tblastn', 'tBLASTn prediction')], max_length=255),
        ),
        migrations.AlterField(
            model_name='fshift',
            name='origin',
            field=models.CharField(choices=[('annotation', 'GenBank annotation'), ('genetack', 'GeneTack prediction'), ('tblastn', 'tBLASTn prediction')], max_length=255),
        ),
    ]
