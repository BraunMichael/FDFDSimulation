function varargout=retreshape(varargin);
% Au del� d'une certaine version de Matlab, le reshape de matrices creuses � plus de 2 dimensions conduit � une erreur
% meme si les dimensions au dela de 2 sont 1
% retreshape transforme alors ces martices en matrices pleines
try;[varargout{1:nargout}]=reshape(varargin{:});
catch;varargin{1}=full(varargin{1});[varargout{1:nargout}]=reshape(varargin{:});end;
