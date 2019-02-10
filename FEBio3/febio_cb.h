#pragma once

class FEModel;

// callback to update window title
bool update_console_cb(FEModel* pfem, unsigned int nwhen, void* pd);

// callback that manages the break points
bool break_point_cb(FEModel* pfem, unsigned int nwhen, void* pd);

// callback for ctrl+c interruptions
bool interrupt_cb(FEModel* pfem, unsigned int nwhen, void* pd);
