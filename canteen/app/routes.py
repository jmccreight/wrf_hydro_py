from flask import render_template, flash, redirect
from app import app
from app.forms import LoginForm, DomainForm
import wrfhydropy


@app.route('/')
@app.route('/index')
def index():
    user = {'username': 'Miguel'}
    posts = [
        {
            'author': {'username': 'John'},
            'body': 'Beautiful day in Portland!'
        },
        {
            'author': {'username': 'Susan'},
            'body': 'The Avengers movie was so cool!'
        }
    ]
    return render_template('index.html', title='Home', user=user, posts=posts)


@app.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        flash('Login requested for user {}, remember_me={}'.format(
            form.username.data, form.remember_me.data))
        return redirect('/index')
    return render_template('login.html', title='Sign In', form=form)


@app.route('/domain', methods=['GET', 'POST'])
def domain():
    form = DomainForm()
    if form.validate_on_submit():
        the_domain = wrfhydropy.Domain(
            domain_top_dir=form.domain_top_dir.data,
            domain_config=form.domain_config.data,
            compatible_version=form.compatible_version.data,
            hydro_namelist_patch_file=form.hydro_namelist_patch_file.data,
            hrldas_namelist_patch_file=form.hrldas_namelist_patch_file.data
        )
        print(the_domain)
        return redirect('/index')
    return render_template('domain.html', title='Init a Domain', form=form)
