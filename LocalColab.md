  FEniCS and Local Runtime in Google Colab
  Contributor: Tamara Dancheva
  
  fenics_tag=2019.1.0.r3
  
  docker run --name notebook-local -w /home/fenics -v $(pwd):/home/fenics/shared -ti -d -p 127.0.0.1:8888:8888 quay.io/fenicsproject/stable:${fenics_tag} "sudo pip install jupyter_http_over_ws; sudo apt-get install -y gmsh; jupyter serverextension enable --py jupyter_http_over_ws; jupyter-notebook --ip=0.0.0.0 --NotebookApp.allow_origin='https://colab.research.google.com' --NotebookApp.port_retries=0 --NotebookApp.allow_root=True --NotebookApp.disable_check_xsrf=True --NotebookApp.token='' --NotebookApp.password='' --port=8888"
  
  FEniCS docker images are available at https://quay.io/repository/fenicsproject/stable?tab=tags
  
  Useful commands:
  
  1. Check the status of notebook-local
  docker logs notebook-local
  If everything is correct, you should see something like below and it is ready to connect
  [W 15:26:49.150 NotebookApp] All authentication is disabled.  Anyone who can connect to this server will be able to run code. 
  jupyter_http_over_ws extension initialized. Listening on /http_over_websocket
  [I 15:26:49.188 NotebookApp] JupyterLab extension loaded from /usr/local/lib/python3.6/dist-packages/jupyterlab
  [I 15:26:49.188 NotebookApp] JupyterLab application directory is /usr/local/share/jupyter/lab
  [I 15:26:49.190 NotebookApp] Serving notebooks from local directory: /home/fenics
  [I 15:26:49.190 NotebookApp] The Jupyter Notebook is running at:
  [I 15:26:49.190 NotebookApp] http://(6d58312f59e0 or 127.0.0.1):8888/
  [I 15:26:49.190 NotebookApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
  [W 15:26:49.195 NotebookApp] No web browser found: could not locate runnable browser.