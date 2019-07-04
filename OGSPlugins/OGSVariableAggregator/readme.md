## OGS Variable Aggregator

The **OGS Variable Aggregator** lets the user obtain composite variables in a similar manner as descrived in the _VarDescriptor_ from _bit.sea_. At the moment, only scalar arrays are accepted in this filter.

The most common variables are already defined. They are:

* _P_c_ = _P1c_ + _P2c_ + _P3c_ + _P4c_
* _P_l_ = _P1l_ + _P2l_ + _P3l_ + _P4l_
* _T_c_ = _B1c_ + _P1c_ + _P2c_ + _P3c_ + _P4c_ + _Z3c_ + _Z4c_ + _Z5c_ + _Z6c_
* _RTc_ = _resMEZ1c_ + _resMEZ2c_ + _resMIZ1c_ + _resMIZ2c_ + _resPBAc_ + _resPPY1c_ + _resPPY2c_ + _resPPY3c_ + _resPPY4c_

The user can generate variables in two ways:

1. By inputting an XML file conforming to the _VarDescriptor_ standard.
2. By using the text box to define the aggregate variable using XML, e.g., 
```xml
<aggregate>
<aggvar name="PP_c">
<var name="P1c"></var> 
<var name="P2c"></var>
</aggvar>
</aggregate>
```

Internally, the XML is parsed using the [pugixml](https://pugixml.org/) library.