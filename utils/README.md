![Screen shot](https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/4299/versions/7/screenshot.jpg)
# timestr
String representation of time without date in MATLAB

Convert serial date numbers into time strings. Time strings are of HH:MM:SS.SSS format - like date strings, but without the date.
Ex:
timestr(now)
09:15:36.4000
Note: A time datetime object was introduced in MATLAB R2014b which is way better than this little utility. If you've got access to it, use it instead!
