extern crate byteorder;

use byteorder::{ReadBytesExt, WriteBytesExt, BigEndian, LittleEndian};
use std::fs::File;
use std::io::prelude::*;

fn main() {
    let mut buffer = File::create("./test_rs.ppm").unwrap();
    // let f: f32 = 2.0;
    // buffer.write(&[f as u8]);
    // wtr.write(b"P6\n");
    // wtr.write_u64::<LittleEndian>(width as u64).unwrap();
    // // wtr.write(width);
    // wtr.write(b" ");
    // wtr.write_u64::<LittleEndian>(height as u64).unwrap();
    // wtr.write(b"\n255\n");
    buffer.write(b"P6\n");
    // let mut wtr = vec![];
    // wtr.write_u64::<LittleEndian>(640).unwrap();
    // wtr.write_u16::<BigEndian>(640).unwrap();
    // buffer.write(&wtr);
    let w: usize = 640;
    // buffer.write(b"{}");
    write!(buffer, "{}", w);
    buffer.write(b" ");
    let h: usize = 480;
    write!(buffer, "{}", h);
    // buffer.write(b"480");
    // let mut wtr = vec![];
    // wtr.write_u64::<LittleEndian>(480).unwrap();
    // wtr.write_u16::<BigEndian>(480).unwrap();
    // buffer.write(&wtr);
    buffer.write(b"\n255\n");
}
