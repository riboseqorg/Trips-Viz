"""
Plugins to add behavior to mpld3 charts
=======================================

Plugins are means of adding additional javascript features to D3-rendered
matplotlib plots.  A number of plugins are defined here; it is also possible
to create nearly any imaginable behavior by defining your own custom plugin.
"""

__all__ = ['connect', 'clear', 'get_plugins', 'PluginBase',
		   'Reset', 'Zoom', 'BoxZoom',
		   'PointLabelTooltip', 'PointHTMLTooltip', 'LineLabelTooltip',
		   'MousePosition']

import collections
import json
import uuid
import matplotlib

from .utils import get_id


def get_plugins(fig):
	"""Get the list of plugins in the figure"""
	connect(fig)
	return fig.mpld3_plugins


def connect(fig, *plugins):
	"""Connect one or more plugins to a figure

	Parameters
	----------
	fig : matplotlib Figure instance
		The figure to which the plugins will be connected

	*plugins :
		Additional arguments should be plugins which will be connected
		to the figure.

	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import plugins
	>>> fig, ax = plt.subplots()
	>>> lines = ax.plot(range(10), '-k')
	>>> plugins.connect(fig, plugins.LineLabelTooltip(lines[0]))
	"""
	if not isinstance(fig, matplotlib.figure.Figure):
		raise ValueError("plugins.connect: first argument must be a figure")
	if not hasattr(fig, 'mpld3_plugins'):
		fig.mpld3_plugins = DEFAULT_PLUGINS[:]
	for plugin in plugins:
		fig.mpld3_plugins.append(plugin)


def clear(fig):
	"""Clear all plugins from the figure, including defaults"""
	fig.mpld3_plugins = []


class PluginBase(object):
	def get_dict(self):
		return self.dict_

	def javascript(self):
		if hasattr(self, "JAVASCRIPT"):
			if hasattr(self, "js_args_"):
				return self.JAVASCRIPT.render(self.js_args_)
			else:
				return self.JAVASCRIPT
		else:
			return ""

	def css(self):
		if hasattr(self, "css_"):
			return self.css_
		else:
			return ""


class Reset(PluginBase):
	"""A Plugin to add a reset button"""
	dict_ = {"type": "reset"}


class MousePosition(PluginBase):
	"""A Plugin to display coordinates for the current mouse position

	Example
	-------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> plugins.connect(fig, plugins.MousePosition())
	>>> fig_to_html(fig)
	"""

	def __init__(self, fontsize=12, fmt=".3g"):
		self.dict_ = {"type": "mouseposition",
					  "fontsize": fontsize,
					  "fmt": fmt}


class Zoom(PluginBase):
	"""A Plugin to add zoom behavior to the plot

	Parameters
	----------
	button : boolean, optional
		if True (default), then add a button to enable/disable zoom behavior
	enabled : boolean, optional
		specify whether the zoom should be enabled by default. By default,
		zoom is enabled if button == False, and disabled if button == True.

	Notes
	-----
	Even if ``enabled`` is specified, other plugins may modify this state.
	"""
	def __init__(self, button=True, enabled=None):
		if enabled is None:
			enabled = not button
		self.dict_ = {"type": "zoom",
					  "button": button,
					  "enabled": enabled}


class BoxZoom(PluginBase):
	"""A Plugin to add box-zoom behavior to the plot

	Parameters
	----------
	button : boolean, optional
		if True (default), then add a button to enable/disable zoom behavior
	enabled : boolean, optional
		specify whether the zoom should be enabled by default. By default,
		zoom is enabled if button == False, and disabled if button == True.

	Notes
	-----
	Even if ``enabled`` is specified, other plugins may modify this state.
	"""
	def __init__(self, button=True, enabled=None):
		if enabled is None:
			enabled = not button
		self.dict_ = {"type": "boxzoom",
					  "button": button,
					  "enabled": enabled}









"""
mpld3.register_plugin("downloadprofile", DownloadProfile);
mpld3_Plugin.prototype.draw = function() {};
mpld3.DownloadProfile = mpld3_DownloadProfile;
mpld3.register_plugin("reset", mpld3_DownloadProfile);
mpld3_DownloadProfile.prototype = Object.create(mpld3_Plugin.prototype);
mpld3_DownloadProfile.prototype.constructor = mpld3_DownloadProfile;
mpld3_DownloadProfile.prototype.requiredProps = [];
mpld3_DownloadProfile.prototype.defaultProps = {};
function mpld3_DownloadProfile(fig, props) {
mpld3_Plugin.call(this, fig, props);
var ResetButton = mpld3.ButtonFactory({
buttonID: "reset",
sticky: false,
onActivate: function() {
	this.toolbar.fig.reset();
},
icon: function() {
	return mpld3.icons["reset"];
}
});
this.fig.buttons.push(ResetButton);
}
"""











'''
	function mpld3_DownloadProfile(fig, props) {
	mpld3_Plugin.call(this, fig, props);
	var ResetButton = mpld3.ButtonFactory({
	buttonID: "reset",
	sticky: false,
	onActivate: function() {
		this.toolbar.fig.reset();
	},
	icon: function() {
		return mpld3.icons["reset"];
	}
	});
	this.fig.buttons.push(ResetButton);
	}

'''









class DownloadProfile(PluginBase):
	"""A Plugin to enable an HTML tooltip:
	formated text which hovers over points.

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		The figure element to apply the tooltip to
	labels : list
		The labels for each point in points, as strings of unescaped HTML.
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
	>>> plugins.connect(fig, PointHTMLTooltip(points[0], labels))
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """




	/* FileSaver.js
	 * A saveAs() FileSaver implementation.
	 * 1.3.8
	 * 2018-03-22 14:03:47
	 *
	 * By Eli Grey, https://eligrey.com
	 * License: MIT
	 *   See https://github.com/eligrey/FileSaver.js/blob/master/LICENSE.md
	 */

	/*global self */
	/*jslint bitwise: true, indent: 4, laxbreak: true, laxcomma: true, smarttabs: true, plusplus: true */

	/*! @source http://purl.eligrey.com/github/FileSaver.js/blob/master/src/FileSaver.js */

	var saveAs = saveAs || (function(view) {
	"use strict";
	// IE <10 is explicitly unsupported
	if (typeof view === "undefined" || typeof navigator !== "undefined" && /MSIE [1-9]\./.test(navigator.userAgent)) {
		return;
	}
	var
	doc = view.document
	// only get URL when necessary in case Blob.js hasn't overridden it yet
	, get_URL = function() {
		return view.URL || view.webkitURL || view;
	}
	, save_link = doc.createElementNS("http://www.w3.org/1999/xhtml", "a")
	, can_use_save_link = "download" in save_link
	, click = function(node) {
		var event = new MouseEvent("click");
		node.dispatchEvent(event);
	}
	, is_safari = /constructor/i.test(view.HTMLElement) || view.safari
	, is_chrome_ios =/CriOS\/[\d]+/.test(navigator.userAgent)
	, setImmediate = view.setImmediate || view.setTimeout
	, throw_outside = function(ex) {
		setImmediate(function() {
			throw ex;
		}, 0);
	}
	, force_saveable_type = "application/octet-stream"
	// the Blob API is fundamentally broken as there is no "downloadfinished" event to subscribe to
	, arbitrary_revoke_timeout = 1000 * 40 // in ms
	, revoke = function(file) {
		var revoker = function() {
			if (typeof file === "string") { // file is an object URL
				get_URL().revokeObjectURL(file);
			} else { // file is a File
				file.remove();
			}
		};
		setTimeout(revoker, arbitrary_revoke_timeout);
	}
	, dispatch = function(filesaver, event_types, event) {
		event_types = [].concat(event_types);
		var i = event_types.length;
		while (i--) {
			var listener = filesaver["on" + event_types[i]];
			if (typeof listener === "function") {
				try {
					listener.call(filesaver, event || filesaver);
				} catch (ex) {
					throw_outside(ex);
				}
			}
		}
	}
	, auto_bom = function(blob) {
		// prepend BOM for UTF-8 XML and text/* types (including HTML)
		// note: your browser will automatically convert UTF-16 U+FEFF to EF BB BF
		if (/^\s*(?:text\/\S*|application\/xml|\S*\/\S*\+xml)\s*;.*charset\s*=\s*utf-8/i.test(blob.type)) {
			return new Blob([String.fromCharCode(0xFEFF), blob], {type: blob.type});
		}
		return blob;
	}
	, FileSaver = function(blob, name, no_auto_bom) {
		if (!no_auto_bom) {
			blob = auto_bom(blob);
		}
		// First try a.download, then web filesystem, then object URLs
		var
		filesaver = this
		, type = blob.type
		, force = type === force_saveable_type
		, object_url
		, dispatch_all = function() {
			dispatch(filesaver, "writestart progress write writeend".split(" "));
		}
		// on any filesys errors revert to saving with object URLs
		, fs_error = function() {
			if ((is_chrome_ios || (force && is_safari)) && view.FileReader) {
				// Safari doesn't allow downloading of blob urls
				var reader = new FileReader();
				reader.onloadend = function() {
					var url = is_chrome_ios ? reader.result : reader.result.replace(/^data:[^;]*;/, 'data:attachment/file;');
					var popup = view.open(url, '_blank');
					if(!popup) view.location.href = url;
							   url=undefined; // release reference before dispatching
							   filesaver.readyState = filesaver.DONE;
					dispatch_all();
				};
				reader.readAsDataURL(blob);
				filesaver.readyState = filesaver.INIT;
				return;
			}
			// don't create more object URLs than needed
			if (!object_url) {
				object_url = get_URL().createObjectURL(blob);
			}
			if (force) {
				view.location.href = object_url;
			} else {
				var opened = view.open(object_url, "_blank");
				if (!opened) {
					// Apple does not allow window.open, see https://developer.apple.com/library/safari/documentation/Tools/Conceptual/SafariExtensionGuide/WorkingwithWindowsandTabs/WorkingwithWindowsandTabs.html
					view.location.href = object_url;
				}
			}
			filesaver.readyState = filesaver.DONE;
			dispatch_all();
			revoke(object_url);
		}
		;
		filesaver.readyState = filesaver.INIT;

		if (can_use_save_link) {
			object_url = get_URL().createObjectURL(blob);
			setImmediate(function() {
				save_link.href = object_url;
				save_link.download = name;
				click(save_link);
				dispatch_all();
				revoke(object_url);
				filesaver.readyState = filesaver.DONE;
			}, 0);
			return;
		}

		fs_error();
	}
	, FS_proto = FileSaver.prototype
	, saveAs = function(blob, name, no_auto_bom) {
		return new FileSaver(blob, name || blob.name || "download", no_auto_bom);
	}
	;

	// IE 10+ (native saveAs)
	if (typeof navigator !== "undefined" && navigator.msSaveOrOpenBlob) {
		return function(blob, name, no_auto_bom) {
			name = name || blob.name || "download";

			if (!no_auto_bom) {
				blob = auto_bom(blob);
			}
			return navigator.msSaveOrOpenBlob(blob, name);
		};
	}

	// todo: detect chrome extensions & packaged apps
	//save_link.target = "_blank";

	FS_proto.abort = function(){};
	FS_proto.readyState = FS_proto.INIT = 0;
	FS_proto.WRITING = 1;
	FS_proto.DONE = 2;

	FS_proto.error =
	FS_proto.onwritestart =
	FS_proto.onprogress =
	FS_proto.onwrite =
	FS_proto.onabort =
	FS_proto.onerror =
	FS_proto.onwriteend =
	null;

	return saveAs;
}(
	typeof self !== "undefined" && self
	|| typeof window !== "undefined" && window
	|| this
));





	var my_icon =
	"data:image/png;base64,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";


	mpld3.register_plugin("downloadprofile", DownloadProfile);
	DownloadProfile.prototype = Object.create(mpld3.Plugin.prototype);

	DownloadProfile.prototype.constructor = DownloadProfile;
	DownloadProfile.prototype.requiredProps = ["returnstr"];
	DownloadProfile.prototype.defaultProps = {};



	function DownloadProfile(fig, props){
		mpld3.Plugin.call(this, fig, props);
		var n = (this.props.returnstr).toString();



		var ResetButton = mpld3.ButtonFactory({
		buttonID: "reset",
		sticky: false,
		onActivate: function() {



			var blob = new Blob([n], {
				type: "text/plain;charset=utf-8;",
			});
			saveAs(blob, "Profile.csv");

		},
		icon: function(){return my_icon;
		}
		});
		this.fig.buttons.push(ResetButton);






	};

	"""

	def __init__(self,returnstr=2017,css=None):
		self.returnstr = returnstr
		self.css_ = css or ""

		self.dict_ = {"type": "downloadprofile",
					  "returnstr": returnstr}





class DownloadPNG(PluginBase):
	"""A Plugin to enable an HTML tooltip:
	formated text which hovers over points.

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		The figure element to apply the tooltip to
	labels : list
		The labels for each point in points, as strings of unescaped HTML.
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
	>>> plugins.connect(fig, PointHTMLTooltip(points[0], labels))
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """




	/* FileSaver.js
	 * A saveAs() FileSaver implementation.
	 * 1.3.8
	 * 2018-03-22 14:03:47
	 *
	 * By Eli Grey, https://eligrey.com
	 * License: MIT
	 *   See https://github.com/eligrey/FileSaver.js/blob/master/LICENSE.md
	 */

	/*global self */
	/*jslint bitwise: true, indent: 4, laxbreak: true, laxcomma: true, smarttabs: true, plusplus: true */

	/*! @source http://purl.eligrey.com/github/FileSaver.js/blob/master/src/FileSaver.js */

	var saveAs = saveAs || (function(view) {
	"use strict";
	// IE <10 is explicitly unsupported
	if (typeof view === "undefined" || typeof navigator !== "undefined" && /MSIE [1-9]\./.test(navigator.userAgent)) {
		return;
	}
	var
	doc = view.document
	// only get URL when necessary in case Blob.js hasn't overridden it yet
	, get_URL = function() {
		return view.URL || view.webkitURL || view;
	}
	, save_link = doc.createElementNS("http://www.w3.org/1999/xhtml", "a")
	, can_use_save_link = "download" in save_link
	, click = function(node) {
		var event = new MouseEvent("click");
		node.dispatchEvent(event);
	}
	, is_safari = /constructor/i.test(view.HTMLElement) || view.safari
	, is_chrome_ios =/CriOS\/[\d]+/.test(navigator.userAgent)
	, setImmediate = view.setImmediate || view.setTimeout
	, throw_outside = function(ex) {
		setImmediate(function() {
			throw ex;
		}, 0);
	}
	, force_saveable_type = "application/octet-stream"
	// the Blob API is fundamentally broken as there is no "downloadfinished" event to subscribe to
	, arbitrary_revoke_timeout = 1000 * 40 // in ms
	, revoke = function(file) {
		var revoker = function() {
			if (typeof file === "string") { // file is an object URL
				get_URL().revokeObjectURL(file);
			} else { // file is a File
				file.remove();
			}
		};
		setTimeout(revoker, arbitrary_revoke_timeout);
	}
	, dispatch = function(filesaver, event_types, event) {
		event_types = [].concat(event_types);
		var i = event_types.length;
		while (i--) {
			var listener = filesaver["on" + event_types[i]];
			if (typeof listener === "function") {
				try {
					listener.call(filesaver, event || filesaver);
				} catch (ex) {
					throw_outside(ex);
				}
			}
		}
	}
	, auto_bom = function(blob) {
		// prepend BOM for UTF-8 XML and text/* types (including HTML)
		// note: your browser will automatically convert UTF-16 U+FEFF to EF BB BF
		if (/^\s*(?:text\/\S*|application\/xml|\S*\/\S*\+xml)\s*;.*charset\s*=\s*utf-8/i.test(blob.type)) {
			return new Blob([String.fromCharCode(0xFEFF), blob], {type: blob.type});
		}
		return blob;
	}
	, FileSaver = function(blob, name, no_auto_bom) {
		if (!no_auto_bom) {
			blob = auto_bom(blob);
		}
		// First try a.download, then web filesystem, then object URLs
		var
		filesaver = this
		, type = blob.type
		, force = type === force_saveable_type
		, object_url
		, dispatch_all = function() {
			dispatch(filesaver, "writestart progress write writeend".split(" "));
		}
		// on any filesys errors revert to saving with object URLs
		, fs_error = function() {
			if ((is_chrome_ios || (force && is_safari)) && view.FileReader) {
				// Safari doesn't allow downloading of blob urls
				var reader = new FileReader();
				reader.onloadend = function() {
					var url = is_chrome_ios ? reader.result : reader.result.replace(/^data:[^;]*;/, 'data:attachment/file;');
					var popup = view.open(url, '_blank');
					if(!popup) view.location.href = url;
							   url=undefined; // release reference before dispatching
							   filesaver.readyState = filesaver.DONE;
					dispatch_all();
				};
				reader.readAsDataURL(blob);
				filesaver.readyState = filesaver.INIT;
				return;
			}
			// don't create more object URLs than needed
			if (!object_url) {
				object_url = get_URL().createObjectURL(blob);
			}
			if (force) {
				view.location.href = object_url;
			} else {
				var opened = view.open(object_url, "_blank");
				if (!opened) {
					// Apple does not allow window.open, see https://developer.apple.com/library/safari/documentation/Tools/Conceptual/SafariExtensionGuide/WorkingwithWindowsandTabs/WorkingwithWindowsandTabs.html
					view.location.href = object_url;
				}
			}
			filesaver.readyState = filesaver.DONE;
			dispatch_all();
			revoke(object_url);
		}
		;
		filesaver.readyState = filesaver.INIT;

		if (can_use_save_link) {
			object_url = get_URL().createObjectURL(blob);
			setImmediate(function() {
				save_link.href = object_url;
				save_link.download = name;
				click(save_link);
				dispatch_all();
				revoke(object_url);
				filesaver.readyState = filesaver.DONE;
			}, 0);
			return;
		}

		fs_error();
	}
	, FS_proto = FileSaver.prototype
	, saveAs = function(blob, name, no_auto_bom) {
		return new FileSaver(blob, name || blob.name || "download", no_auto_bom);
	}
	;

	// IE 10+ (native saveAs)
	if (typeof navigator !== "undefined" && navigator.msSaveOrOpenBlob) {
		return function(blob, name, no_auto_bom) {
			name = name || blob.name || "download";

			if (!no_auto_bom) {
				blob = auto_bom(blob);
			}
			return navigator.msSaveOrOpenBlob(blob, name);
		};
	}

	// todo: detect chrome extensions & packaged apps
	//save_link.target = "_blank";

	FS_proto.abort = function(){};
	FS_proto.readyState = FS_proto.INIT = 0;
	FS_proto.WRITING = 1;
	FS_proto.DONE = 2;

	FS_proto.error =
	FS_proto.onwritestart =
	FS_proto.onprogress =
	FS_proto.onwrite =
	FS_proto.onabort =
	FS_proto.onerror =
	FS_proto.onwriteend =
	null;

	return saveAs;
}(
	typeof self !== "undefined" && self
	|| typeof window !== "undefined" && window
	|| this
));





	var my_save_icon = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAADIAAAAyCAYAAAAeP4ixAAAABHNCSVQICAgIfAhkiAAAAitJREFUaIHtms+rqkAcxY+PSpAgcOVKgtzopqD2LfuTw12b9gWDbWYzuLFNIIjSD+gtLleuVtZXm95V32enzMz3HDgzo6PK9Xq9ogH8+dcC3kVjjHTyN4QQ8DwPcRzf7WDbNhzHKV1wuVzeHXs6nWI4HAIAPM/DbrcDAHS7XczncwwGg8JxM0aEEFitVoUdTqcTQfYtcRwjiqKb+5fLJVPjZxvXdbFYLArNZKLleV4lkbJIkgSu6yIMw4dtMkYexek38MxMrSZ7kZlaGQEem6mdEeDLTH5RqqURADifz5nr2hrJc7Mhysa27bt7ka7rlcb9uJEqTwVFtDdah8MBnHMZWjI1qJCN+L4P3/fJhWTTmGj9N/LbaIyRlya7qqrodD6+5QD4euE6Ho9P272kbjabwbKsyqLKwDnHer1+2q5d0crzfRSmKAp+HospivIeVSUgG9lsNthutzK0pIzHY0wmE1KfxkSrMUbI0bJtG6PRSIaWlF6vR+5DNqKqKlRVJReSTXujxTmX/j5iWRZ5AyYbiaII+/2e2o2EYRjkPu2Nlmma6Pf7MrSklDlRIRvRdb3y0Y0M2hutIAgQBIEMLSmGYZAnfCkjsh8aAfrK1d5o6bou/VnrI6uWaZowTZNcSDbtilYYhtJXqqLar5Axomna3W/gjDEwxt6j7E1ompa5zkRL1rcLGeS1Kvm/g4QQYIwhSZKPCnsVTdPgOE76u8c3N0bqSmNWrcYY+Qu5IscJiMF9TQAAAABJRU5ErkJggg==";


	mpld3.register_plugin("downloadpng", DownloadPNG);
	DownloadPNG.prototype = Object.create(mpld3.Plugin.prototype);

	DownloadPNG.prototype.constructor = DownloadPNG;
	DownloadPNG.prototype.requiredProps = ["returnstr"];
	DownloadPNG.prototype.defaultProps = {};


	function saveCanvas(x_canvas){
		x_canvas.toBlob(function(blob) {
			saveAs(
				blob
				, "screenshot.png"
			);
		}, "image/png");
	}


	function DownloadPNG(fig, props){
		mpld3.Plugin.call(this, fig, props);
		var n = (this.props.returnstr).toString();



		var ResetButton = mpld3.ButtonFactory({
		buttonID: "reset",
		sticky: false,
		onActivate: function() {



			//var blob = new Blob([n], {
			//	type: "text/plain;charset=utf-8;",
			//});
			//saveAs(blob, "Image.html");

			//base_image = new Image();
			//base_image.src = 'data:image/png;base64,'+n
			//var canvas = document.createElement('canvas');
			//canvas.id     = "YourCanvas";
			//var canvas = document.getElementById('YourCanvas');
			//context = canvas.getContext('2d');
			// Draw image within
			//context.drawImage(base_image, 0,0);
			//canvas.toBlob(function(blob) {
			//	saveAs(blob, "Image.png");
			//}, "image/png");

			base_image = new Image();
			base_image.src ='data:image/png;base64,'+n;
			base_image.onload = function(){
				var canvas = document.createElement('canvas');
				canvas.width = 2300;
				canvas.height = 1100;
				context = canvas.getContext('2d');

				// Draw image within
				context.drawImage(base_image, 0,0);
				// Save the canvas
				saveCanvas(canvas);
				// Remove it
			};





		},
		icon: function(){return my_save_icon;
		}
		});
		this.fig.buttons.push(ResetButton);






	};

	"""

	def __init__(self,returnstr=2017,css=None):
		self.returnstr = returnstr
		self.css_ = css or ""

		self.dict_ = {"type": "downloadpng",
					  "returnstr": returnstr}










class PointLabelTooltip(PluginBase):
	"""A Plugin to enable a tooltip: text which hovers over points.

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		The figure element to apply the tooltip to
	labels : array or None
		If supplied, specify the labels for each point in points.  If not
		supplied, the (x, y) values will be used.
	hoffset, voffset : integer
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10

	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> plugins.connect(fig, PointLabelTooltip(points[0]))
	>>> fig_to_html(fig)
	"""
	def __init__(self, points, labels=None,
				 hoffset=0, voffset=10, location="mouse"):
		if location not in ["bottom left", "top left", "bottom right",
							"top right", "mouse"]:
			raise ValueError("invalid location: {0}".format(location))
		if isinstance(points, matplotlib.lines.Line2D):
			suffix = "pts"
		else:
			suffix = None
		self.dict_ = {"type": "tooltip",
					  "id": get_id(points, suffix),
					  "labels": labels,
					  "hoffset": hoffset,
					  "voffset": voffset,
					  "location": location}


class LineLabelTooltip(PluginBase):
	"""A Plugin to enable a tooltip: text which hovers over a line.

	Parameters
	----------
	line : matplotlib Line2D object
		The figure element to apply the tooltip to
	label : string
		If supplied, specify the labels for each point in points.  If not
		supplied, the (x, y) values will be used.
	hoffset, voffset : integer
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10

	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> lines = ax.plot(range(10), 'o')
	>>> plugins.connect(fig, LineLabelTooltip(lines[0]))
	>>> fig_to_html(fig)
	"""
	def __init__(self, points, label=None,
				 hoffset=0, voffset=10, location="mouse"):
		if location not in ["bottom left", "top left", "bottom right",
							"top right", "mouse"]:
			raise ValueError("invalid location: {0}".format(location))
		self.dict_ = {"type": "tooltip",
					  "id": get_id(points),
					  "labels": label if label is None else [label],
					  "hoffset": hoffset,
					  "voffset": voffset,
					  "location": location}


class LinkedBrush(PluginBase):
	"""A Plugin to enable linked brushing between plots

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		A representative of the scatter plot elements to brush.
	button : boolean, optional
		if True (default), then add a button to enable/disable zoom behavior
	enabled : boolean, optional
		specify whether the zoom should be enabled by default. default=True.

	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> import numpy as np
	>>> from mpld3 import fig_to_html, plugins
	>>> X = np.random.random((3, 100))
	>>> fig, ax = plt.subplots(3, 3)
	>>> for i in range(2):
	...     for j in range(2):
	...         points = ax[i, j].scatter(X[i], X[j])
	>>> plugins.connect(fig, LinkedBrush(points))
	>>> fig_to_html(fig)

	Notes
	-----
	Notice that in the above example, only one of the four sets of points is
	passed to the plugin. This is all that is needed: for the sake of efficient
	data storage, mpld3 keeps track of which plot objects draw from the same
	data.

	Also note that for the linked brushing to work correctly, the data must
	not contain any NaNs. The presence of NaNs makes the different data views
	have different sizes, so that mpld3 is unable to link the related points.
	"""

	def __init__(self, points, button=True, enabled=True):
		if isinstance(points, matplotlib.lines.Line2D):
			suffix = "pts"
		else:
			suffix = None

		self.dict_ = {"type": "linkedbrush",
					  "button": button,
					  "enabled": enabled,
					  "id": get_id(points, suffix)}


class PointHTMLTooltip(PluginBase):
	"""A Plugin to enable an HTML tooltip:
	formated text which hovers over points.

	Parameters
	----------
	points : matplotlib Collection or Line2D object
		The figure element to apply the tooltip to
	labels : list
		The labels for each point in points, as strings of unescaped HTML.
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> labels = ['<h1>{title}</h1>'.format(title=i) for i in range(10)]
	>>> plugins.connect(fig, PointHTMLTooltip(points[0], labels))
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """
	mpld3.register_plugin("htmltooltip", HtmlTooltipPlugin);
	HtmlTooltipPlugin.prototype = Object.create(mpld3.Plugin.prototype);
	HtmlTooltipPlugin.prototype.constructor = HtmlTooltipPlugin;
	HtmlTooltipPlugin.prototype.requiredProps = ["id"];
	HtmlTooltipPlugin.prototype.defaultProps = {labels:null,
												hoffset:0,
												voffset:10};
	function HtmlTooltipPlugin(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	HtmlTooltipPlugin.prototype.draw = function(){
	   var obj = mpld3.get_element(this.props.id);
	   var labels = this.props.labels;
	   var tooltip = d3.select("body").append("div")
					.attr("class", "mpld3-tooltip")
					.style("position", "absolute")
					.style("z-index", "10")
					.style("visibility", "hidden");

	   obj.elements()
		   .on("mouseover", function(d, i){
							  tooltip.html(labels[i])
									 .style("visibility", "visible");})
		   .on("mousemove", function(d, i){
				  tooltip
					.style("top", d3.event.pageY + this.props.voffset + "px")
					.style("left",d3.event.pageX + this.props.hoffset + "px");
				 }.bind(this))
		   .on("mouseout",  function(d, i){
						   tooltip.style("visibility", "hidden");});
	};
	"""

	def __init__(self, points, labels=None,
				 hoffset=0, voffset=10, css=None):
		self.points = points
		self.labels = labels
		self.voffset = voffset
		self.hoffset = hoffset
		self.css_ = css or ""
		if isinstance(points, matplotlib.lines.Line2D):
			suffix = "pts"
		else:
			suffix = None
		self.dict_ = {"type": "htmltooltip",
					  "id": get_id(points, suffix),
					  "labels": labels,
					  "hoffset": hoffset,
					  "voffset": voffset}


class LineHTMLTooltip(PluginBase):
	"""A Plugin to enable an HTML tooltip:
	formated text which hovers over points.

	Parameters
	----------
	points : matplotlib Line2D object
		The figure element to apply the tooltip to
	label : string
		The label for the line, as strings of unescaped HTML.
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> lines = ax.plot(range(10))
	>>> label = '<h1>line {title}</h1>'.format(title='A')
	>>> plugins.connect(fig, LineHTMLTooltip(lines[0], label))
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """
	mpld3.register_plugin("linehtmltooltip", LineHTMLTooltip);
	LineHTMLTooltip.prototype = Object.create(mpld3.Plugin.prototype);
	LineHTMLTooltip.prototype.constructor = LineHTMLTooltip;
	LineHTMLTooltip.prototype.requiredProps = ["id"];
	LineHTMLTooltip.prototype.defaultProps = {label:null,
											  hoffset:0,
											  voffset:10};
	function LineHTMLTooltip(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	LineHTMLTooltip.prototype.draw = function(){
		var obj = mpld3.get_element(this.props.id, this.fig);
		var label = this.props.label
		var tooltip = d3.select("body").append("div")
					.attr("class", "mpld3-tooltip")
					.style("position", "absolute")
					.style("z-index", "10")
					.style("visibility", "hidden");

		obj.elements()
		   .on("mouseover", function(d, i){
							   tooltip.html(label)
									  .style("visibility", "visible");
									 })
			.on("mousemove", function(d, i){
				  tooltip
					.style("top", d3.event.pageY + this.props.voffset + "px")
					.style("left",d3.event.pageX + this.props.hoffset + "px");
				 }.bind(this))
		   .on("mouseout",  function(d, i){
						   tooltip.style("visibility", "hidden");})
	};
	"""

	def __init__(self, line, label=None,
				 hoffset=0, voffset=10,
				 css=None):
		self.line = line
		self.label = label
		self.voffset = voffset
		self.hoffset = hoffset
		self.css_ = css or ""
		self.dict_ = {"type": "linehtmltooltip",
					  "id": get_id(line),
					  "label": label,
					  "hoffset": hoffset,
					  "voffset": voffset}









class InteractiveLegendPlugin(PluginBase):
	"""A plugin for an interactive legends.

	Inspired by http://bl.ocks.org/simzou/6439398

	Parameters
	----------
	plot_elements : iterable of matplotliblib elements
		the elements to associate with a given legend items
	labels : iterable of strings
		The labels for each legend element
	ax :  matplotlib axes instance, optional
		the ax to which the legend belongs. Default is the first
		axes. The legend will be plotted to the right of the specified
		axes
	alpha_sel : float, optional
		the alpha value to apply to the plot_element(s) associated
		with the legend item when the legend item is selected.
		Default is 1.0
	alpha_unsel : float, optional
		the alpha value to apply to the plot_element(s) associated
		with the legend item when the legend item is unselected.
		Default is 0.2
	xoffset : int, optional
		apply x offset to the rectangles and labels
	yoffset : int, optional
		apply y offset to the rectangles and labels
	start_visible : boolean, optional (could be a list of booleans)
		defines if objects should start selected on not.
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld4 import fig_to_html, plugins
	>>> N_paths = 5
	>>> N_steps = 100
	>>> x = np.linspace(0, 10, 100)
	>>> y = 0.1 * (np.random.random((N_paths, N_steps)) - 0.5)
	>>> y = y.cumsum(1)
	>>> fig, ax = plt.subplots()
	>>> labels = ["a", "b", "c", "d", "e"]
	>>> line_collections = ax.plot(x, y.T, lw=4, alpha=0.1)
	>>> interactive_legend = plugins.InteractiveLegendPlugin(line_collections,
	...                                                      labels,
	...                                                      alpha_unsel=0.1)
	>>> plugins.connect(fig, interactive_legend)
	>>> fig_to_html(fig)
	"""

	JAVASCRIPT = """
	mpld3.register_plugin("interactive_legend", InteractiveLegend);
	InteractiveLegend.prototype = Object.create(mpld3.Plugin.prototype);
	InteractiveLegend.prototype.constructor = InteractiveLegend;
	InteractiveLegend.prototype.requiredProps = ["element_ids", "labels"];
	InteractiveLegend.prototype.defaultProps = {"ax":null,
												"alpha_sel":1.0,
												"alpha_unsel":0,
												"xoffset":0,
												"yoffset":0,
												"start_visible":[false,false,false,false,false,false,false,false,false,false,false]}
	function InteractiveLegend(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	InteractiveLegend.prototype.draw = function(){
		console.log(this);
		var alpha_sel = this.props.alpha_sel;
		var alpha_unsel = this.props.alpha_unsel;
		var xoffset = this.props.xoffset;
		var yoffset = this.props.yoffset;

		var legendItems = new Array();
		for(var i=0; i<this.props.labels.length; i++){
			var obj = {};
			obj.label = this.props.labels[i];

			var element_id = this.props.element_ids[i];
			mpld3_elements = [];
			for(var j=0; j<element_id.length; j++){
				var mpld3_element = mpld3.get_element(element_id[j], this.fig);

				// mpld3_element might be null in case of Line2D instances
				// for we pass the id for both the line and the markers. Either
				// one might not exist on the D3 side
				if(mpld3_element){
					mpld3_elements.push(mpld3_element);
				}
			}

			obj.mpld3_elements = mpld3_elements;
			console.log("HHHHHHHEEEEEEEELLLLLLLLPPPPPPPP"+i+this.props.start_visible[i])
			obj.visible = this.props.start_visible[i]; // should become be setable from python side
			legendItems.push(obj);

		}
		console.log("LEGEND ITEMS:"+legendItems);

		// determine the axes with which this legend is associated
		var ax = this.props.ax
		if(!ax){
			ax = this.fig.axes[0];
		} else{
			ax = mpld3.get_element(ax, this.fig);
		}

		// add a legend group to the canvas of the figure
		var legend = this.fig.canvas.append("svg:g")
							   .attr("class", "legend");

		// add the cds_areagles
		legend.selectAll("rect")
				.data(legendItems)
			 .enter().append("rect")
				.attr("height",10)
				.attr("width", 20)
				.attr("x",ax.width+1+ax.position[0]-(xoffset))
				.attr("y",function(d,i) {
							return ax.position[1]+ i * 25 + (yoffset);}) // should be -10
				.attr("stroke", get_color)
				.attr("class", "legend-box")
				.style("fill", function(d, i) {
							return d.visible ? get_color(d) : "white";})
				.on("click", click);

		// add the labels
		legend.selectAll("text")
			  .data(legendItems)
			  .enter().append("text")
			  .attr("x", function (d) {
							return ax.width+5+ax.position[0] - (xoffset-25);}) // should be +25
			  .attr("y", function(d,i) {
							return ax.position[1]+ (i * 25)+(yoffset+10);})
			  .attr("fontsize",1)
			  .text(function(d) { return d.label });

		// specify the action on click
		function click(d,i){
			d.visible = !d.visible;
			d3.select(this)
			  .style("fill",function(d, i) {
				return d.visible ? get_color(d) : "white";
			  })

			for(var i=0; i<d.mpld3_elements.length; i++){
				var type = d.mpld3_elements[i].constructor.name;
				if(type =="mpld3_Line"){
					d3.select(d.mpld3_elements[i].path[0][0])
						.style("stroke-opacity",
								d.visible ? alpha_sel : alpha_unsel);
				} else if((type=="mpld3_PathCollection")||
						 (type=="mpld3_Markers")){
					d3.selectAll(d.mpld3_elements[i].pathsobj[0])
						.style("stroke-opacity",
								d.visible ? alpha_sel : alpha_unsel)
						.style("fill-opacity",
								d.visible ? alpha_sel : alpha_unsel);
				} else{
					console.log(type + " not yet supported");
				}
			}
		};

		// helper function for determining the color of the cds_areagles
		function get_color(d){
			var type = d.mpld3_elements[0].constructor.name;
			var color = "black";
			if(type =="mpld3_Line"){
				color = d.mpld3_elements[0].props.edgecolor;
			} else if((type=="mpld3_PathCollection")||
					  (type=="mpld3_Markers")){
				color = d.mpld3_elements[0].props.facecolors[0];
			} else{
				console.log(type + " not yet supported");
			}
			return color;
		};
	};
	"""

	css_ = """
	.legend-box {
	  cursor: pointer;
	}
	"""

	def __init__(self, plot_elements, labels, ax=None,
				 alpha_sel=1, alpha_unsel=0.2, xoffset=0, yoffset=0,start_visible=False):

		self.ax = ax

		if ax:
			ax = get_id(ax)

		# start_visible could be a list
		if isinstance(start_visible, bool):
			start_visible = [start_visible] * len(labels)
		elif not len(start_visible) == len(labels):
			raise ValueError("{} out of {} visible params has been set".format(len(start_visible), len(labels)))

		mpld4_element_ids = self._determine_mpld4ids(plot_elements)
		self.mpld4_element_ids = mpld4_element_ids
		self.dict_ = {"type": "interactive_legend",
					  "element_ids": mpld4_element_ids,
					  "labels": labels,
					  "ax": ax,
					  "alpha_sel": alpha_sel,
					  "alpha_unsel": alpha_unsel,
					  "xoffset":xoffset,
					  "yoffset":yoffset,
					  "start_visible":start_visible}

	def _determine_mpld4ids(self, plot_elements):
		"""
		Helper function to get the mpld4_id for each
		of the specified elements.
		"""
		mpld4_element_ids = []

		# There are two things being done here. First,
		# we make sure that we have a list of lists, where
		# each inner list is associated with a single legend
		# item. Second, in case of Line2D object we pass
		# the id for both the marker and the line.
		# on the javascript side we filter out the nulls in
		# case either the line or the marker has no equivalent
		# D3 representation.
		for entry in plot_elements:
			ids = []
			if isinstance(entry, collections.Iterable):
				for element in entry:
					mpld4_id = get_id(element)
					ids.append(mpld4_id)
					if isinstance(element, matplotlib.lines.Line2D):
						mpld4_id = get_id(element, 'pts')
						ids.append(mpld4_id)
			else:
				ids.append(get_id(entry))
				if isinstance(entry, matplotlib.lines.Line2D):
					mpld4_id = get_id(entry, 'pts')
					ids.append(mpld4_id)
			mpld4_element_ids.append(ids)
		return mpld4_element_ids


class PointClickableHTMLTooltip(PluginBase):
	"""A plugin for pop-up windows with data with rich HTML

	Parameters
	----------
	points : matplotlib Collection object
		The figure element to apply the tooltip to
	labels : list
		The labels for each point in points, as strings of unescaped HTML.
	targets : list
		The target data or rich HTML to be displayed when each collection element is clicked
	hoffset, voffset : integer, optional
		The number of pixels to offset the tooltip text.  Default is
		hoffset = 0, voffset = 10
	css : str, optional
		css to be included, for styling the label html and target data/tables, if desired
	Examples
	--------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import plugins
	>>> fig, ax = plt.subplots(1,1)
	>>> xx = yy = range(10)
	>>> scat = ax.scatter(xx, range(10))
	>>> targets = map(lambda (x, y): "<marquee>It works!<br><h1>{}, {}</h1></marquee>".format(x, y),
	>>>               zip(xx, yy))
	>>> labels = map(lambda (x, y): "{}, {}".format(x,y), zip(xx, yy))
	>>> from mpld3.plugins import PointClickableHTMLTooltip
	>>> plugins.connect(fig, PointClickableHTMLTooltip(scat, labels=labels, targets=targets))

	"""

	JAVASCRIPT="""
	mpld3.register_plugin("clickablehtmltooltip", PointClickableHTMLTooltip);
	PointClickableHTMLTooltip.prototype = Object.create(mpld3.Plugin.prototype);
	PointClickableHTMLTooltip.prototype.constructor = PointClickableHTMLTooltip;
	PointClickableHTMLTooltip.prototype.requiredProps = ["id"];
	PointClickableHTMLTooltip.prototype.defaultProps = {labels:null,
												 targets:null,
												 hoffset:0,
												 voffset:10};
	function PointClickableHTMLTooltip(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	PointClickableHTMLTooltip.prototype.draw = function(){
	   var obj = mpld3.get_element(this.props.id);
	   var labels = this.props.labels;
	   var targets = this.props.targets;

	   var tooltip = d3.select("body").append("div")
					.attr("class", "mpld3-tooltip")
					.style("position", "absolute")
					.style("z-index", "10")
					.style("visibility", "hidden");

	   obj.elements()
		   .on("mouseover", function(d, i){
				  if ($(obj.elements()[0][0]).css( "fill-opacity" ) > 0.5 || $(obj.elements()[0][0]).css( "stroke-opacity" ) > 0) {
							  tooltip.html(labels[i])
									 .style("visibility", "visible");
							  } })

		   .on("mousedown", function(d, i){
							  window.open().document.write(targets[i]);
							   })
		   .on("mousemove", function(d, i){
				  tooltip
					.style("top", d3.event.pageY + this.props.voffset + "px")
					.style("left",d3.event.pageX + this.props.hoffset + "px");
				 }.bind(this))
		   .on("mouseout",  function(d, i){
						   tooltip.style("visibility", "hidden");});
	};
	"""
	def __init__(self, points, labels=None, targets=None,
				 hoffset=2, voffset=-6, css=None):
		self.points = points
		self.labels = labels
		self.targets = targets
		self.voffset = voffset
		self.hoffset = hoffset
		self.css_ = css or ""
		if targets is not None:
			styled_targets = map(lambda x: self.css_ + x, targets)
		else:
			styled_targets = None


		if isinstance(points, matplotlib.lines.Line2D):
			suffix = "pts"
		else:
			suffix = None
		self.dict_ = {"type": "clickablehtmltooltip",
					  "id": get_id(points, suffix),
					  "labels": labels,
					  "targets": styled_targets,
					  "hoffset": hoffset,
					  "voffset": voffset}


class MouseXPosition(PluginBase):
	"""Like MousePosition, but only show the X coordinate"""

	JAVASCRIPT="""
  mpld3.register_plugin("mousexposition", MouseXPositionPlugin);
  MouseXPositionPlugin.prototype = Object.create(mpld3.Plugin.prototype);
  MouseXPositionPlugin.prototype.constructor = MouseXPositionPlugin;
  MouseXPositionPlugin.prototype.requiredProps = [];
  MouseXPositionPlugin.prototype.defaultProps = {
	fontsize: 12,
	fmt: "0d"
  };
  function MouseXPositionPlugin(fig, props) {
	mpld3.Plugin.call(this, fig, props);
  }
  MouseXPositionPlugin.prototype.draw = function() {
	var fig = this.fig;
	var fmt = d3.format(this.props.fmt);
	var coords = fig.canvas.append("text").attr("class", "mpld3-coordinates").style("text-anchor", "end").style("font-size", this.props.fontsize).attr("x", this.fig.width - 5).attr("y", this.fig.height - 5);
	for (var i = 0; i < this.fig.axes.length; i++) {
	  var update_coords = function() {
		var ax = fig.axes[i];
		return function() {
		  var pos = d3.mouse(this), x = ax.x.invert(pos[0]), y = ax.y.invert(pos[1]);
		  coords.text(fmt(x));
		};
	  }();
	  fig.axes[i].baseaxes.on("mousemove", update_coords).on("mouseout", function() {
		coords.text("");
	  });
	}
  };"""
	"""A Plugin to display coordinates for the current mouse position

	Example
	-------
	>>> import matplotlib.pyplot as plt
	>>> from mpld3 import fig_to_html, plugins
	>>> fig, ax = plt.subplots()
	>>> points = ax.plot(range(10), 'o')
	>>> plugins.connect(fig, plugins.MouseXPosition())
	>>> fig_to_html(fig)
	"""

	def __init__(self, fontsize=12, fmt="8.0f"):
		self.dict_ = {"type": "mousexposition",
					  "fontsize": fontsize,
					  "fmt": fmt}


class TopToolbar(PluginBase):
	"""Plugin for moving toolbar to top of figure"""

	JAVASCRIPT = """
	mpld3.register_plugin("toptoolbar", TopToolbar);
	TopToolbar.prototype = Object.create(mpld3.Plugin.prototype);
	TopToolbar.prototype.constructor = TopToolbar;
	TopToolbar.prototype.defaultProps = {"xoffset":0,
										 "yoffset":0}
	function TopToolbar(fig, props){
		mpld3.Plugin.call(this, fig, props);
	};

	TopToolbar.prototype.draw = function(){
	  // the toolbar svg doesn't exist
	  // yet, so first draw it
	  this.fig.toolbar.draw();
	  var xoffset = this.props.xoffset;
	  var yoffset = this.props.yoffset;

	  // set the x and y position of the toolbar
	  this.fig.toolbar.toolbar.attr("y", 800+yoffset);
	  this.fig.toolbar.toolbar.attr("x", 910+xoffset);

	  // set the width and height of the toolbar
	  this.fig.toolbar.toolbar.attr("width", 280);
	  this.fig.toolbar.toolbar.attr("height", 30);

	  // set the width and height of the buttons
	  this.fig.toolbar.buttonsobj.attr("width",24);
	  this.fig.toolbar.buttonsobj.attr("height",24);

	  // set the space between the buttons
	  this.fig.toolbar.buttonsobj.attr("x", function(d, i) {return i * 30;});

	  // show toolbar
	  this.fig.toolbar.buttonsobj.transition(0).attr("y", 0);

	  // remove event triggers so that toolbar doesnt dissapear when mouse leaves plot.
	  this.fig.canvas
		.on("mouseenter", null)
		.on("mouseleave", null)
		.on("touchenter", null)
		.on("touchstart", null);

	  // then remove the draw function,
	  // so that it is not called again
	  this.fig.toolbar.draw = function() {}
	}
	"""
	def __init__(self,xoffset=0, yoffset=0):
		self.dict_ = {"type": "toptoolbar",
					  "xoffset":xoffset,
					  "yoffset":yoffset}


DEFAULT_PLUGINS = [Reset(), Zoom(), BoxZoom()]
