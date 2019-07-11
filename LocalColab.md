  # Local Runtimes in Google Colab
  Contributor: Tamara Dancheva and Van-Dang Nguyen
  
  Google Colab notebooks can be connected to either a hosted runtime provided by Google Cloud or a local runtime. The hosted runtime allows us to access free resources for up to 12 hours at a time. For longer excecutions, it is more convenient to connect to the local runtimes. We consider two options: direct local runtime and local runtime with FEniCS containers.

![#f03c15](https://placehold.it/15/f03c15/000000?text=+) `WARNINGS`

(1) It is very important to know that this option allows you to execute code on your local hardware and have access to your local file system. Please be careful when removing folders or files.

(2) We are not reponsible for any data loss due to the use of this functionality.
 
  ## Direct local runtime
  This option allows for a direct connection to a local computer. The full instruction is available here https://research.google.com/colaboratory/local-runtimes.html
  ```bash
  jupyter notebook --NotebookApp.allow_origin='https://colab.research.google.com' --port=8886 --NotebookApp.port_retries=0
  ```
  ## Local runtime with FEniCS containers
  This option allows for connecting to existing FEniCS containers provided by https://fenicsproject.org.
  ```bash
  fenics_tag=2019.1.0.r3
  docker run --name notebook-local -w /home/fenics -v $(pwd):/home/fenics/shared -ti -d -p 127.0.0.1:8888:8888 quay.io/fenicsproject/stable:${fenics_tag} "sudo pip install jupyter_http_over_ws; sudo apt-get install -y gmsh; jupyter serverextension enable --py jupyter_http_over_ws; jupyter-notebook --ip=0.0.0.0 --NotebookApp.allow_origin='https://colab.research.google.com' --NotebookApp.port_retries=0 --NotebookApp.allow_root=True --NotebookApp.disable_check_xsrf=True --NotebookApp.token='' --NotebookApp.password='' --port=8888"
  ```
  FEniCS docker tags are available at https://quay.io/repository/fenicsproject/stable?tab=tags
  
  ### Useful commands:
  
  ##### Check the status of notebook-local
  ```bash
  docker logs notebook-local
  ```
  If everything is correct, you should see something like below
  ```bash
  [W 15:26:49.150 NotebookApp] All authentication is disabled.  Anyone who can connect to this server will be able to run code. 
  jupyter_http_over_ws extension initialized. Listening on /http_over_websocket
  [I 15:26:49.188 NotebookApp] JupyterLab extension loaded from /usr/local/lib/python3.6/dist-packages/jupyterlab
  [I 15:26:49.188 NotebookApp] JupyterLab application directory is /usr/local/share/jupyter/lab
  [I 15:26:49.190 NotebookApp] Serving notebooks from local directory: /home/fenics
  [I 15:26:49.190 NotebookApp] The Jupyter Notebook is running at:
  [I 15:26:49.190 NotebookApp] http://(6d58312f59e0 or 127.0.0.1):8888/
  [I 15:26:49.190 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
  [W 15:26:49.195 NotebookApp] No web browser found: could not locate runnable browser.
  ```
 It is now ready to connect: Google Colab > CONNECT > Connect to local runtime...> Backend port 8888 > CONNECT.
 
 ##### List Docker containers
  ```bash
  docker ps
  ```
 ##### Stop and remove an existing containers
  ```bash
  docker stop <CONTAINER ID>
  docker rm <CONTAINER ID>
  ```
