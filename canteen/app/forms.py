from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, BooleanField, SubmitField
from wtforms.validators import DataRequired

class LoginForm(FlaskForm):
    username = StringField('Username', validators=[DataRequired()])
    password = PasswordField('Password', validators=[DataRequired()])
    remember_me = BooleanField('Remember Me')
    submit = SubmitField('Sign In')


class DomainForm(FlaskForm):
    domain_top_dir = StringField(
        'domain_top_dir',
        validators=[DataRequired()]
    )
    domain_config = StringField(
        'domain_config',
        validators=[DataRequired()]
    )
    compatible_version = StringField(
        'compatible_version',
        validators=[DataRequired()],
        default='None'
    )
    hydro_namelist_patch_file = StringField(
        'hydro_namelist_patch_file',
        validators=[DataRequired()],
        default='hydro_namelist_patches.json'
    )
    hrldas_namelist_patch_file = StringField(
        'hrldas_namelist_patch_file',
        validators=[DataRequired()],
        default='hrldas_namelist_patches.json'
    )
    submit = SubmitField('Initalize')
