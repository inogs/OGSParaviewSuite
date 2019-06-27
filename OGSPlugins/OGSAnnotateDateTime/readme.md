## OGS Annotate Date Time

The **OGS Annotate Date Time** plugin is a simple tool to print in the GUI the current instant date and time.

The **OGS Reader** uses the time structure to store the time stamp in the time array. By default, this plugin reads the data from the timestamp and converts it to a readable format using the _ctime_ library. At the user's request, the metadata array can be used instead to obtain a representation of the time. This solution is usually more accurate since it involves less conversions. The output is in the form the user has inputed using the [date time format](http://strftime.org/).

Some of the most used formats are included below as reference:

| Code | Meaning           | Example  |
| ---- |------------------ | -------- |
| %d   | Day of the month  | 30       |
| %m   | Month of the year | 12       |
| %b   | Month of the year | Dec      |
| %B   | Month of the year | December |
| %y   | Year              | 19       |
| %Y   | Year              | 2019     |
| %H   | Hour (24h clock)  | 20       |
| %I   | Hour (12h clock)  | 08       |
| %M   | Minute            | 15       |
| %S   | Second            | 00       |
| %p   | AM or PM          | AM       |
