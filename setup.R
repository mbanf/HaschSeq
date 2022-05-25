# install.packages("renv")
# renv::init()

# setup R v3.3.2.
# Install it
$ git clone https://github.com/jcrodriguez1989/renv-installer.git ~/.renv
$ echo 'export RENV_ROOT="$HOME/.renv"' >> ~/.bashrc
$ echo 'export PATH="$RENV_ROOT/bin:$PATH"' >> ~/.bashrc
$ echo -e 'if command -v renv 1>/dev/null 2>&1; then\n  eval "$(renv init -)"\nfi' >> ~/.bashrc
$ exec "$SHELL"

# Install user-locally R 3.0.0
$ renv install 3.4
# Set it to use it when R is called in this folder
$ renv local 3.3.2
$ R